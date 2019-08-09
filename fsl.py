"""
Analysis workflows (FSL) for ds000109.

This script replicates, using three distinct software pipelines, the contrast shown in
figure 5 of https://www.jneurosci.org/content/32/16/5553 ; which compare two task conditions
over two age groups (young / old and false-belief / false-photo)

"""

# General Set up ###
from nipype.pipeline import engine as pe
from nipype.algorithms.modelgen import SpecifyModel
from nipype.interfaces import fsl, utility as niu, io as nio
from nipype.workflows.fmri.fsl.preprocess import create_susan_smooth
from niworkflows.interfaces.bids import DerivativesDataSink as BIDSDerivatives
from interfaces import PtoZ
from .utils import _dict_ds, _bids2nipypeinfo, _len, _dof, _neg


DATA_ITEMS = ['bold', 'mask', 'events', 'regressors', 'tr']


class DerivativesDataSink(BIDSDerivatives):
    out_path_base = 'FSLAnalysis'


class GroupDerivativesDataSink(BIDSDerivatives):
    out_path_base = 'FSL-all'


def first_level_wf(in_files, output_dir, fwhm=6.0, name='fsl_1st_level'):
    """
    Creates the first level of analysis (individual participants,
    concatenating runs).

    We aim to reproduce 2 contrast maps. The first is for "young" subjects,
    and the second "old" subjects.
    Both are the same contrast, false_belief_story vs. false_belief_photo
    """
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

    # Preprocessing fMRIPrep is reluctant to do: SUSAN smoothing
    susan = create_susan_smooth()
    susan.inputs.inputnode.fwhm = fwhm

    # l1_spec adapts the "Info" dictionary from BIDS2info to the
    # FSL tooling.
    l1_spec = pe.Node(SpecifyModel(
        parameter_source='FSL',
        input_units='secs',
        high_pass_filter_cutoff=100
    ), name='l1_spec')

    # l1_model creates a first-level model design
    l1_model = pe.Node(fsl.Level1Design(
        bases={'dgamma': {'derivs': True}},
        model_serial_correlations=True,
        contrasts=[('story>photo', 'T',
                    ['false_belief_story', 'false_belief_photo'],
                    [1, -1])],
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
        desc='storyvsphoto'), name='ds_cope', run_without_submitting=True)

    ds_varcope = pe.Node(DerivativesDataSink(
        base_directory=str(output_dir), keep_dtype=False, suffix='varcope',
        desc='storyvsphoto'), name='ds_varcope', run_without_submitting=True)

    ds_zstat = pe.Node(DerivativesDataSink(
        base_directory=str(output_dir), keep_dtype=False, suffix='zstat',
        desc='storyvsphoto'), name='ds_zstat', run_without_submitting=True)

    ds_tstat = pe.Node(DerivativesDataSink(
        base_directory=str(output_dir), keep_dtype=False, suffix='tstat',
        desc='storyvsphoto'), name='ds_tstat', run_without_submitting=True)

    # Connect up the workflow
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


def second_level_wf(output_dir, bids_ref, name='wf_2nd_level'):
    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['group_mask', 'in_copes', 'in_varcopes']),
        name='inputnode')

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['zstats_raw', 'zstats_fwe', 'zstats_clust',
                'clust_index_file', 'clust_localmax_txt_file']),
        name='outputnode')

    # Configure FSL 2nd level analysis
    l2_model = pe.Node(fsl.L2Model(), name='l2_model')
    flameo = pe.Node(fsl.FLAMEO(run_mode='flame12'), name='flameo')

    merge_copes = pe.Node(fsl.Merge(dimension='t'), name='merge_copes')
    merge_varcopes = pe.Node(fsl.Merge(dimension='t'), name='merge_varcopes')

    # Thresholding - FWE ################################################
    # smoothest -r %s -d %i -m %s
    smoothness = pe.Node(fsl.SmoothEstimate(), name='smoothness')
    # ptoz 0.025 -g %f
    # p = 0.05 / 2 for 2-tailed test
    fwe_ptoz = pe.Node(PtoZ(pvalue=0.025), name='fwe_ptoz')
    # fslmaths %s -uthr %s -thr %s nonsignificant
    # fslmaths %s -sub nonsignificant zstat1_thresh
    fwe_nonsig0 = pe.Node(fsl.Threshold(direction='above'), name='fwe_nonsig0')
    fwe_nonsig1 = pe.Node(fsl.Threshold(direction='below'), name='fwe_nonsig1')
    fwe_thresh = pe.Node(fsl.BinaryMaths(operation='sub'), name='fwe_thresh')

    # Thresholding - Cluster ############################################
    # cluster -i %s -c %s -t 3.2 -p 0.025 -d %s --volume=%s  \
    #     --othresh=thresh_cluster_fwe_zstat1 --connectivity=26 --mm
    cluster_kwargs = {
        'connectivity': 26,
        'threshold': 3.2,
        'pthreshold': 0.025,
        'out_threshold_file': True,
        'out_index_file': True,
        'out_localmax_txt_file': True
    }
    cluster_pos = pe.Node(fsl.Cluster(**cluster_kwargs),
                          name='cluster_pos')
    cluster_neg = pe.Node(fsl.Cluster(**cluster_kwargs),
                          name='cluster_neg')
    zstat_inv = pe.Node(fsl.BinaryMaths(operation='mul', operand_value=-1),
                        name='zstat_inv')
    cluster_inv = pe.Node(fsl.BinaryMaths(operation='mul', operand_value=-1),
                          name='cluster_inv')
    cluster_all = pe.Node(fsl.BinaryMaths(operation='add'), name='cluster_all')

    ds_zraw = pe.Node(GroupDerivativesDataSink(
        base_directory=str(output_dir), keep_dtype=False, suffix='zstat', sub='all'),
        name='ds_zraw', run_without_submitting=True)
    ds_zraw.inputs.source_file = bids_ref

    ds_zfwe = pe.Node(GroupDerivativesDataSink(
        base_directory=str(output_dir), keep_dtype=False, suffix='zstat',
        desc='fwe', sub='all'), name='ds_zfwe', run_without_submitting=True)
    ds_zfwe.inputs.source_file = bids_ref

    ds_zclust = pe.Node(GroupDerivativesDataSink(
        base_directory=str(output_dir), keep_dtype=False, suffix='zstat',
        desc='clust', sub='all'), name='ds_zclust', run_without_submitting=True)
    ds_zclust.inputs.source_file = bids_ref

    ds_clustidx_pos = pe.Node(GroupDerivativesDataSink(
        base_directory=str(output_dir), keep_dtype=False, suffix='pclusterindex', sub='all'),
        name='ds_clustidx_pos', run_without_submitting=True)
    ds_clustidx_pos.inputs.source_file = bids_ref

    ds_clustlmax_pos = pe.Node(GroupDerivativesDataSink(
        base_directory=str(output_dir), keep_dtype=False, suffix='plocalmax',
        desc='storyvsphoto', sub='all'), name='ds_clustlmax_pos', run_without_submitting=True)
    ds_clustlmax_pos.inputs.source_file = bids_ref

    ds_clustidx_neg = pe.Node(GroupDerivativesDataSink(
        base_directory=str(output_dir), keep_dtype=False, suffix='nclusterindex', sub='all'),
        name='ds_clustidx_neg', run_without_submitting=True)
    ds_clustidx_neg.inputs.source_file = bids_ref

    ds_clustlmax_neg = pe.Node(GroupDerivativesDataSink(
        base_directory=str(output_dir), keep_dtype=False, suffix='nlocalmax',
        desc='storyvsphoto', sub='all'), name='ds_clustlmax_neg', run_without_submitting=True)
    ds_clustlmax_neg.inputs.source_file = bids_ref

    workflow.connect([
        (inputnode, l2_model, [(('in_copes', _len), 'num_copes')]),
        (inputnode, flameo, [('group_mask', 'mask_file')]),
        (inputnode, smoothness, [('group_mask', 'mask_file'),
                                 (('in_copes', _dof), 'dof')]),
        (inputnode, merge_copes, [('in_copes', 'in_files')]),
        (inputnode, merge_varcopes, [('in_varcopes', 'in_files')]),

        (l2_model, flameo, [('design_mat', 'design_file'),
                            ('design_con', 't_con_file'),
                            ('design_grp', 'cov_split_file')]),
        (merge_copes, flameo, [('merged_file', 'cope_file')]),
        (merge_varcopes, flameo, [('merged_file', 'var_cope_file')]),
        (flameo, smoothness, [('res4d', 'residual_fit_file')]),

        (flameo, fwe_nonsig0, [('zstats', 'in_file')]),
        (fwe_nonsig0, fwe_nonsig1, [('out_file', 'in_file')]),
        (smoothness, fwe_ptoz, [('resels', 'resels')]),
        (fwe_ptoz, fwe_nonsig0, [('zstat', 'thresh')]),
        (fwe_ptoz, fwe_nonsig1, [(('zstat', _neg), 'thresh')]),
        (flameo, fwe_thresh, [('zstats', 'in_file')]),
        (fwe_nonsig1, fwe_thresh, [('out_file', 'operand_file')]),

        (flameo, cluster_pos, [('zstats', 'in_file')]),
        (merge_copes, cluster_pos, [('merged_file', 'cope_file')]),
        (smoothness, cluster_pos, [('volume', 'volume'),
                                   ('dlh', 'dlh')]),
        (flameo, zstat_inv, [('zstats', 'in_file')]),
        (zstat_inv, cluster_neg, [('out_file', 'in_file')]),
        (cluster_neg, cluster_inv, [('threshold_file', 'in_file')]),
        (merge_copes, cluster_neg, [('merged_file', 'cope_file')]),
        (smoothness, cluster_neg, [('volume', 'volume'),
                                   ('dlh', 'dlh')]),
        (cluster_pos, cluster_all, [('threshold_file', 'in_file')]),
        (cluster_inv, cluster_all, [('out_file', 'operand_file')]),

        (flameo, ds_zraw, [('zstats', 'in_file')]),
        (fwe_thresh, ds_zfwe, [('out_file', 'in_file')]),
        (cluster_all, ds_zclust, [('out_file', 'in_file')]),
        (cluster_pos, ds_clustidx_pos, [('index_file', 'in_file')]),
        (cluster_pos, ds_clustlmax_pos, [('localmax_txt_file', 'in_file')]),
        (cluster_neg, ds_clustidx_neg, [('index_file', 'in_file')]),
        (cluster_neg, ds_clustlmax_neg, [('localmax_txt_file', 'in_file')]),
    ])
    return workflow
