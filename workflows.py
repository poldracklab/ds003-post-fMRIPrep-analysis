"""
Analysis workflows
"""

from nipype.pipeline import engine as pe
from nipype.algorithms.modelgen import SpecifyModel
from nipype.interfaces import fsl, utility as niu, io as nio
from nipype.workflows.fmri.fsl.preprocess import create_susan_smooth
from niworkflows.interfaces.bids import DerivativesDataSink as BIDSDerivatives


class DerivativesDataSink(BIDSDerivatives):
    out_path_base = 'FSLAnalysis'


def first_level_wf(output_dir, fwhm=6.0, name='wf_1st_level'):
    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['in_file', 'in_mask', 'events_file', 'regressors', 'repetition_time']),
        name='inputnode')

    # Extract motion parameters from regressors file
    runinfo = pe.Node(niu.Function(
        input_names=['in_file', 'events_file', 'regressors_file', 'regressors_names'],
        function=_bids2nipypeinfo, output_names=['info', 'realign_file']),
        name='runinfo')

    # Set the column names to be used from the confounds file
    runinfo.inputs.regressors_names = ['global_signal', 'dvars', 'framewise_displacement'] + \
        ['a_comp_cor_%02d' % i for i in range(6)]

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
    feat_fit = pe.Node(fsl.FEAT(), name='feat_fit')

    feat_select = pe.Node(nio.SelectFiles({
        'cope': 'cope1.nii.gz',
        'pe': 'pe[0-9][0-9].nii.gz',
        'tstat': 'tstat1.nii.gz',
        'varcope': 'varcope1.nii.gz',
        'zstat': 'zstat1.nii.gz',
    }), name='feat_select')

    ds_cope = pe.Node(DerivativesDataSink(
        base_directory=output_dir, keep_dtype=False, stat='cope', suffix='statmap',
        cont='intask'), name='ds_cope', run_without_submitting=True)

    ds_varcope = pe.Node(DerivativesDataSink(
        base_directory=output_dir, keep_dtype=False, stat='varcope', suffix='statmap',
        cont='intask'), name='ds_varcope', run_without_submitting=True)

    ds_zstat = pe.Node(DerivativesDataSink(
        base_directory=output_dir, keep_dtype=False, stat='z', suffix='statmap',
        cont='intask'), name='ds_zstat', run_without_submitting=True)

    ds_tstat = pe.Node(DerivativesDataSink(
        base_directory=output_dir, keep_dtype=False, stat='t', suffix='statmap',
        cont='intask'), name='ds_tstat', run_without_submitting=True)

    workflow.connect([
        (inputnode, susan, [('in_file', 'inputnode.in_files'),
                            ('in_mask', 'inputnode.mask_file')]),
        (inputnode, runinfo, [
            ('events_file', 'events_file'),
            ('regressors', 'regressors_file')]),
        (susan, l1_spec, [('outputnode.smoothed_files', 'functional_runs')]),
        (inputnode, l1_spec, [('repetition_time', 'time_repetition')]),
        (inputnode, l1_model, [('repetition_time', 'interscan_interval')]),
        (inputnode, ds_cope, [('repetition_time', 'source_file')]),
        (inputnode, ds_varcope, [('repetition_time', 'source_file')]),
        (inputnode, ds_zstat, [('repetition_time', 'source_file')]),
        (inputnode, ds_tstat, [('repetition_time', 'source_file')]),
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
        (feat_select, ds_cope, [('varcope', 'in_file')]),
        (feat_select, ds_cope, [('zstat', 'in_file')]),
        (feat_select, ds_cope, [('tstat', 'in_file')]),
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
