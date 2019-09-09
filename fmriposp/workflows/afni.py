"""
Analysis workflows (AFNI) for ds000109.

This script replicates, using three distinct software pipelines, the contrast shown in
figure 5 of https://www.jneurosci.org/content/32/16/5553 ; which compare two task conditions
over two age groups (young / old and false-belief / false-photo)

"""

# General Set up ###
from nipype.pipeline import engine as pe
from nipype.interfaces import afni
from niworkflows.interfaces.bids import DerivativesDataSink as BIDSDerivatives

DATA_ITEMS = ['bold', 'mask', 'events', 'regressors', 'tr']


class DerivativesDataSink(BIDSDerivatives):
    out_path_base = 'AFNIAnalysis'


class GroupDerivativesDataSink(BIDSDerivatives):
    out_path_base = 'AFNI-all'


def first_level_wf(in_files, output_dir, fwhm=6.0, name='afni_1st_level'):
    """
    Creates the first level of analysis (individual runs).

    We aim to reproduce 2 contrast maps. The first is for "young" subjects,
    and the second "old" subjects.
    Both are the same contrast, false_belief_story vs. false_belief_photo
    """
    workflow = pe.Workflow(name=name)

    datasource = pe.Node(niu.Function(function=_dict_ds, output_names=DATA_ITEMS),
                         name='datasource')
    datasource.inputs.in_dict = in_files
    datasource.iterables = ('sub', sorted(in_files.keys()))

    # Preprocessing fMRIPrep is reluctant to do: Smoothing using AFNI's 3dmerge tool. 
    # This is the same smoothing tool preferred by the afni_proc.py which AFNI 
    # recommends that analysts use to prepare their pipelines.
    smoothing = pe.Node(afni.Merge(
    	blurfwhm=fwhm
    ), name='smoothing')

    # Use 3dDeconvolve to set up the design matrix file required by 3dRemlFit
    # The only output we need is the x1D matrix. The x1D_stop command halts
    # 3dDeconvolve as soon as the x1D matrix is generated.

    # As inputs to 3dDeconvolve we still need to provide:
    #	 a mask via mask='PATH' (I have a line in the .connect that does this, I think?)

    #	 an indication of how many stimulus files we have via num_stimts=X
    # 	 stimulus timeseries information via a file and stim_times='PATH'
    #	 stimulus labels via a file and stim_label='PATH'
    modelspec = pe.Node(afni.deconvolve(
    	x1D_stop=True,
    ), name='modelspec')

    # Remlfit is used to generate subject-level t-statistics
    # Remlfit needs a 3D+time data set and the x matrix from 3dDeconvolve
    lv1_afni = pe.Node(afni.remlfit(
    	outputtype='NIFTI',
    ), name='lv1_afni')

    # Connect the nodes into the workflow
    workflow.connect([
    	(datasource, smoothing, [('bold', 'input_files')]),

    	#feed smoothed 3d+t data into model spec
    	(smoothing, modelspec, [('out_file', 'in_files')]),

    	#feed smoothed 3d+t data and design matrix into lv1 analysis
    	(smoothing, lv1_afni, [('out_file', 'in_files')]),
    	(modelspec, lv1_afni, [('x1D', 'matrix')]),

    	#If we have a brainmask from Prep we should feed it into various steps:
    	(datasource, modelspec, [('...', 'mask')]),
    	(datasource, lv1_afni, [('...', 'mask')]),

    	#files I think we need to sink
    	(lv1_afni, ..., [
    		('out_file', ...),
    		('wherr_file', ...),])
    ])
    return workflow