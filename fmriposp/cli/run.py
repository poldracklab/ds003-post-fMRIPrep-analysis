#!/usr/bin/env python3
import sys
import logging
import json
from pathlib import Path
from templateflow.api import get as tpl_get, templates as get_tpl_list


logging.addLevelName(25, 'IMPORTANT')  # Add a new level between INFO and WARNING
logging.addLevelName(15, 'VERBOSE')  # Add a new level between INFO and DEBUG
logger = logging.getLogger('workflow')


metadata = {
    'Name': 'ds003 example postprocessing',
    'BIDSVersion': '1.1.1',
    'PipelineDescription': {
        'Name': 'ds003-post-fMRIPrep-analysis'
    },
    'CodeURL': 'https://github.com/poldracklab/ds003-post-fMRIPrep-analysis'
}


def get_parser():
    """Define the command line interface"""
    from argparse import ArgumentParser
    from argparse import RawTextHelpFormatter
    from .. import __version__

    parser = ArgumentParser(description='DS000003 Analysis Workflow',
                            formatter_class=RawTextHelpFormatter)

    # Arguments as specified by BIDS-Apps
    # required, positional arguments
    # IMPORTANT: they must go directly with the parser object
    parser.add_argument(
        'derivatives_dir', action='store', type=Path,
        help='the root folder of a derivatives set generated with fMRIPrep '
             '(sub-XXXXX folders should be found at the top level in this folder).')
    parser.add_argument('output_dir', action='store', type=Path,
                        help='the output path for the outcomes of preprocessing and visual '
                             'reports')
    parser.add_argument('analysis_level', choices=['participant', 'group'], nargs='+',
                        help='processing stage to be run, "participant" means individual analysis '
                             'and "group" is second level analysis.')

    parser.add_argument('--version', action='version', version=__version__)

    # Options that affect how pyBIDS is configured
    g_bids = parser.add_argument_group('Options for filtering BIDS queries')
    g_bids.add_argument('--participant-label', action='store', type=str,
                        nargs='*', help='process only particular subjects')
    g_bids.add_argument('--task', action='store', type=str, nargs='*',
                        help='select a specific task to be processed')
    g_bids.add_argument('--run', action='store', type=int, nargs='*',
                        help='select a specific run identifier to be processed')
    g_bids.add_argument('--space', action='store', choices=get_tpl_list() + ['T1w', 'template'],
                        help='select a specific space to be processed')
    g_bids.add_argument('--bids-dir', action='store', type=Path,
                        help='point to the BIDS root of the dataset from which the derivatives '
                             'were calculated (in case the derivatives folder is not the default '
                             '(i.e. ``BIDS_root/derivatives``).')

    g_perfm = parser.add_argument_group('Options to handle performance')
    g_perfm.add_argument("-v", "--verbose", dest="verbose_count", action="count", default=0,
                         help="increases log verbosity for each occurence, debug level is -vvv")
    g_perfm.add_argument('--ncpus', '--nprocs', action='store', type=int,
                         help='maximum number of threads across all processes')
    g_perfm.add_argument('--nthreads', '--omp-nthreads', action='store', type=int,
                         help='maximum number of threads per-process')

    g_other = parser.add_argument_group('Other options')
    g_other.add_argument('-w', '--work-dir', action='store', type=Path,
                         help='path where intermediate results should be stored')

    return parser


def main():
    """Entry point"""
    from os import cpu_count
    from multiprocessing import set_start_method
    from bids.layout import BIDSLayout
    from nipype import logging as nlogging
    set_start_method('forkserver')

    opts = get_parser().parse_args()

    # Retrieve logging level
    log_level = int(max(25 - 5 * opts.verbose_count, logging.DEBUG))
    # Set logging
    logger.setLevel(log_level)
    nlogging.getLogger('nipype.workflow').setLevel(log_level)
    nlogging.getLogger('nipype.interface').setLevel(log_level)
    nlogging.getLogger('nipype.utils').setLevel(log_level)

    # Resource management options
    plugin_settings = {
        'plugin': 'MultiProc',
        'plugin_args': {
            'n_procs': opts.ncpus,
            ''
            'raise_insufficient': False,
            'maxtasksperchild': 1,
        }
    }
    # Permit overriding plugin config with specific CLI options
    if not opts.ncpus or opts.ncpus < 1:
        plugin_settings['plugin_args']['n_procs'] = cpu_count()

    nthreads = opts.nthreads
    if not nthreads or nthreads < 1:
        nthreads = cpu_count()

    derivatives_dir = opts.derivatives_dir.resolve()
    bids_dir = opts.bids_dir or derivatives_dir.parent

    # Get absolute path to BIDS directory
    bids_dir = opts.bids_dir.resolve()
    layout = BIDSLayout(str(bids_dir), validate=False,
                        ignore=['.git', '.gitattributes', '.datalad'],
                        derivatives=str(derivatives_dir))
    query = {'domains': 'derivatives', 'desc': 'preproc',
             'suffix': 'bold', 'extension': ['.nii', '.nii.gz']}

    if opts.participant_label:
        query['subject'] = opts.participant_label
    if opts.run:
        query['run'] = opts.run
    if opts.task:
        query['task'] = opts.task
    if opts.space:
        query['space'] = opts.space

    # Preprocessed files that are input to the workflow
    prepped_bold = layout.get(**query)
    if not prepped_bold:
        print('No preprocessed files found under the given derivatives '
              'folder "%s".' % derivatives_dir, file=sys.stderr)

    # The magic happens here
    if 'participant' in opts.analysis_level:
        from ..workflows.fsl import participant_level_wf

        output_dir = opts.output_dir.resolve()
        output_dir.mkdir(exist_ok=True, parents=True)
        logger.info('Writting 1st level outputs to "%s".', output_dir)
        base_entities = set(['subject', 'session', 'task', 'run', 'acquisition', 'reconstruction'])
        inputs = {}
        for part in prepped_bold:
            entities = part.entities
            sub = entities['subject']
            logger.info('Gathering preprocessed data for subject %s', sub)
            inputs[sub] = {}
            base = base_entities.intersection(entities)
            subquery = {k: v for k, v in entities.items() if k in base}
            inputs[sub]['bold'] = part.path
            inputs[sub]['mask'] = layout.get(
                domains='derivatives',
                suffix='mask',
                return_type='file',
                extension=['.nii', '.nii.gz'],
                space=query['space'],
                **subquery)
            inputs[sub]['events'] = layout.get(
                suffix='events', return_type='file', **subquery)[0]
            inputs[sub]['regressors'] = layout.get(
                domains='derivatives',
                suffix='regressors',
                return_type='file',
                extension=['.tsv'],
                **subquery)
            inputs[sub]['tr'] = part.get_metadata()['RepetitionTime']

        workflow = participant_level_wf(inputs, output_dir)
        workflow.base_dir = opts.work_dir
        workflow.run(**plugin_settings)

    if 'group' in opts.analysis_level:
        from ..workflows.fsl import group_level_wf
        import re

        output_dir = opts.output_dir.resolve()
        metafile = '{}/FSLAnalysis/dataset_description.json'.format(output_dir)
        with open(metafile, 'w') as metafile:
            json.dump(metadata, metafile, indent=4)
        glayout = BIDSLayout(str(bids_dir), validate=False, derivatives=str(output_dir))

        base_entities = set(['subject', 'session', 'task', 'run', 'acquisition', 'reconstruction'])
        in_copes = []
        in_varcopes = []
        for part in prepped_bold:
            entities = part.entities
            base = base_entities.intersection(entities)
            subquery = {k: v for k, v in entities.items() if k in base}
            in_copes.append(glayout.get(
                domains='derivatives',
                suffix='cope',
                return_type='file',
                extension=['.nii', '.nii.gz'],
                space=query['space'],
                **subquery)[0])
            in_varcopes.append(glayout.get(
                domains='derivatives',
                suffix='varcope',
                return_type='file',
                extension=['.nii', '.nii.gz'],
                space=query['space'],
                **subquery)[0])
        bids_ref = re.sub('sub-[0-9]+', 'sub-all', prepped_bold[0].path)

        group_mask = tpl_get(entities['space'],
                             resolution=2,
                             desc='brain',
                             suffix='mask')

        group_out = output_dir / 'FSLAnalysis'
        group_out.mkdir(exist_ok=True, parents=True)

        workflow = group_level_wf(group_out, bids_ref)

        # set inputs
        workflow.inputs.inputnode.group_mask = str(group_mask)
        workflow.inputs.inputnode.in_copes = in_copes
        workflow.inputs.inputnode.in_varcopes = in_varcopes

        workflow.base_dir = opts.work_dir
        workflow.run(**plugin_settings)

    return 0


if __name__ == '__main__':
    sys.exit(main())
