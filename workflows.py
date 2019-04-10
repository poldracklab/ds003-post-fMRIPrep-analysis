"""
Analysis workflows
"""

from nipype.pipeline import engine as pe
from nipype.algorithms.modelgen import SpecifyModel
from nipype.interfaces import fsl, utility as niu, io as nio
from nipype.workflows.fmri.fsl.preprocess import create_susan_smooth
from niworkflows.interfaces.bids import DerivativesDataSink as BIDSDerivatives


DATA_ITEMS = ['bold', 'mask', 'events', 'regressors', 'tr']


class DerivativesDataSink(BIDSDerivatives):
    out_path_base = 'FSLAnalysis'


def first_level_wf(in_files, output_dir, fwhm=6.0, name='wf_1st_level'):
    workflow = pe.Workflow(name=name)

    datasource = pe.Node(niu.Function(function=_dict_ds, output_names=DATA_ITEMS),
                         name='datasource')
    datasource.inputs.in_dict = in_files
    datasource.iterables = ('sub', sorted(in_files.keys()))

    # Extract motion parameters from regressors file
    runinfo = pe.Node(niu.Function(
        input_names=['in_file', 'events_file', 'regressors_file', 'regressors_names'],
        function=_bids2nipypeinfo, output_names=['info', 'realign_file']),
        name='runinfo')

    # Set the column names to be used from the confounds file
    runinfo.inputs.regressors_names = ['dvars', 'framewise_displacement'] + \
        ['a_comp_cor_%02d' % i for i in range(6)] + ['cosine%02d' % i for i in range(4)]

    # SUSAN smoothing
    susan = create_susan_smooth()
    susan.inputs.inputnode.fwhm = fwhm

    l1_spec = pe.Node(SpecifyModel(
        parameter_source='FSL',
        input_units='secs',
        high_pass_filter_cutoff=100
    ), name='l1_spec')

    # l1_model creates a first-level model design
    l1_model = pe.Node(fsl.Level1Design(
        bases={'dgamma': {'derivs': True}},
        model_serial_correlations=True,
        contrasts=[('intask', 'T', ['word', 'pseudoword'], [1, 1])],
        # orthogonalization=orthogonality,
    ), name='l1_model')

    # feat_spec generates an fsf model specification file
    feat_spec = pe.Node(fsl.FEATModel(), name='feat_spec')
    # feat_fit actually runs FEAT
    feat_fit = pe.Node(fsl.FEAT(), name='feat_fit', mem_gb=12)

    feat_select = pe.Node(nio.SelectFiles({
        'cope': 'stats/cope1.nii.gz',
        'pe': 'stats/pe[0-9][0-9].nii.gz',
        'tstat': 'stats/tstat1.nii.gz',
        'varcope': 'stats/varcope1.nii.gz',
        'zstat': 'stats/zstat1.nii.gz',
    }), name='feat_select')

    ds_cope = pe.Node(DerivativesDataSink(
        base_directory=str(output_dir), keep_dtype=False, suffix='cope',
        desc='intask'), name='ds_cope', run_without_submitting=True)

    ds_varcope = pe.Node(DerivativesDataSink(
        base_directory=str(output_dir), keep_dtype=False, suffix='varcope',
        desc='intask'), name='ds_varcope', run_without_submitting=True)

    ds_zstat = pe.Node(DerivativesDataSink(
        base_directory=str(output_dir), keep_dtype=False, suffix='zstat',
        desc='intask'), name='ds_zstat', run_without_submitting=True)

    ds_tstat = pe.Node(DerivativesDataSink(
        base_directory=str(output_dir), keep_dtype=False, suffix='tstat',
        desc='intask'), name='ds_tstat', run_without_submitting=True)

    workflow.connect([
        (datasource, susan, [('bold', 'inputnode.in_files'),
                             ('mask', 'inputnode.mask_file')]),
        (datasource, runinfo, [
            ('events', 'events_file'),
            ('regressors', 'regressors_file')]),
        (susan, l1_spec, [('outputnode.smoothed_files', 'functional_runs')]),
        (datasource, l1_spec, [('tr', 'time_repetition')]),
        (datasource, l1_model, [('tr', 'interscan_interval')]),
        (datasource, ds_cope, [('bold', 'source_file')]),
        (datasource, ds_varcope, [('bold', 'source_file')]),
        (datasource, ds_zstat, [('bold', 'source_file')]),
        (datasource, ds_tstat, [('bold', 'source_file')]),
        (susan, runinfo, [('outputnode.smoothed_files', 'in_file')]),
        (runinfo, l1_spec, [
            ('info', 'subject_info'),
            ('realign_file', 'realignment_parameters')]),
        (l1_spec, l1_model, [('session_info', 'session_info')]),
        (l1_model, feat_spec, [
            ('fsf_files', 'fsf_file'),
            ('ev_files', 'ev_files')]),
        (l1_model, feat_fit, [('fsf_files', 'fsf_file')]),
        (feat_fit, feat_select, [('feat_dir', 'base_directory')]),
        (feat_select, ds_cope, [('cope', 'in_file')]),
        (feat_select, ds_varcope, [('varcope', 'in_file')]),
        (feat_select, ds_zstat, [('zstat', 'in_file')]),
        (feat_select, ds_tstat, [('tstat', 'in_file')]),
    ])
    return workflow


def second_level_wf(output_dir, name='wf_2nd_level'):
    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['group_mask', 'in_copes', 'in_varcopes']),
        name='inputnode')

    # Configure FSL 2nd level analysis
    l2_model = pe.Node(fsl.L2Model(), name='l2_model')
    flameo_ols = pe.Node(fsl.FLAMEO(run_mode='ols'), name='flameo_ols')

    # Thresholding - FDR ################################################
    # Calculate pvalues with ztop
    fdr_ztop = pe.Node(fsl.ImageMaths(op_string='-ztop', suffix='_pval'),
                       name='fdr_ztop')
    # Find FDR threshold: fdr -i zstat1_pval -m <group_mask> -q 0.05
    # fdr_th = <write Nipype interface for fdr>
    # Apply threshold:
    # fslmaths zstat1_pval -mul -1 -add 1 -thr <fdr_th> -mas <group_mask> \
    #     zstat1_thresh_vox_fdr_pstat1

    # Thresholding - FWE ################################################
    # smoothest -r %s -d %i -m %s
    # ptoz 0.05 -g %f
    # fslmaths %s -thr %s zstat1_thresh

    # Thresholding - Cluster ############################################
    # cluster -i %s -c %s -t 3.2 -p 0.05 -d %s --volume=%s  \
    #     --othresh=thresh_cluster_fwe_zstat1 --connectivity=26 --mm

    workflow.connect([
        (inputnode, l2_model, [(('in_copes', _len), 'num_copes')]),
        (inputnode, flameo_ols, [('group_mask', 'mask_file')]),
        (l2_model, flameo_ols, [('design_mat', 'design_file'),
                                ('design_con', 't_con_file'),
                                ('design_grp', 'cov_split_file')]),
    ])
    return workflow


def _bids2nipypeinfo(in_file, events_file, regressors_file,
                     regressors_names=None,
                     motion_columns=None,
                     decimals=3, amplitude=1.0):
    from pathlib import Path
    import numpy as np
    import pandas as pd
    from nipype.interfaces.base.support import Bunch

    # Process the events file
    events = pd.read_csv(events_file, sep=r'\s+')

    bunch_fields = ['onsets', 'durations', 'amplitudes']

    if not motion_columns:
        from itertools import product
        motion_columns = ['_'.join(v) for v in product(('trans', 'rot'), 'xyz')]

    out_motion = Path('motion.par').resolve()

    regress_data = pd.read_csv(regressors_file, sep=r'\s+')
    np.savetxt(out_motion, regress_data[motion_columns].values, '%g')
    if regressors_names is None:
        regressors_names = sorted(set(regress_data.columns) - set(motion_columns))

    if regressors_names:
        bunch_fields += ['regressor_names']
        bunch_fields += ['regressors']

    runinfo = Bunch(
        scans=in_file,
        conditions=list(set(events.trial_type.values)),
        **{k: [] for k in bunch_fields})

    for condition in runinfo.conditions:
        event = events[events.trial_type.str.match(condition)]

        runinfo.onsets.append(np.round(event.onset.values, 3).tolist())
        runinfo.durations.append(np.round(event.duration.values, 3).tolist())
        if 'amplitudes' in events.columns:
            runinfo.amplitudes.append(np.round(event.amplitudes.values, 3).tolist())
        else:
            runinfo.amplitudes.append([amplitude] * len(event))

    if 'regressor_names' in bunch_fields:
        runinfo.regressor_names = regressors_names
        runinfo.regressors = regress_data[regressors_names].fillna(0.0).values.T.tolist()

    return [runinfo], str(out_motion)


def _get_tr(in_dict):
    return in_dict.get('RepetitionTime')


def _len(inlist):
    return len(inlist)


def _dict_ds(in_dict, sub, order=['bold', 'mask', 'events', 'regressors', 'tr']):
    return tuple([in_dict[sub][k] for k in order])
