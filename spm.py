"""
Analysis workflows (SPM) for ds000109.

This script replicates, using three distinct software pipelines, the contrast shown in
figure 5 of https://www.jneurosci.org/content/32/16/5553 ; which compare two task conditions
over two age groups (young / old and false-belief / false-photo)

"""

# General Set up ###
from nipype.pipeline import engine as pe

DATA_ITEMS = ['bold', 'mask', 'events', 'regressors', 'tr']


class DerivativesDataSink(BIDSDerivatives):
    out_path_base = 'SPMAnalysis'


class GroupDerivativesDataSink(BIDSDerivatives):
    out_path_base = 'SPM-all'


def first_level_wf(in_files, output_dir, fwhm=6.0, name='spm_1st_level'):
    """
    Creates the first level of analysis (individual runs).

    We aim to reproduce 2 contrast maps. The first is for "young" subjects,
    and the second "old" subjects.
    Both are the same contrast, false_belief_story vs. false_belief_photo
    """
    workflow = pe.Workflow(name=name)
