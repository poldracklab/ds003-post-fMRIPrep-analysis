#!/usr/bin/env python3
import sys
import logging
from pathlib import Path
from templateflow.api import templates as get_tpl_list

__version__ = '1.0.0'
logging.addLevelName(25, 'IMPORTANT')  # Add a new level between INFO and WARNING
logging.addLevelName(15, 'VERBOSE')  # Add a new level between INFO and DEBUG
logger = logging.getLogger('cli')


def get_parser():
    """Define the command line interface"""
    from argparse import ArgumentParser
    from argparse import RawTextHelpFormatter

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
    layout = BIDSLayout(str(bids_dir), validate=False, derivatives=str(derivatives_dir))
    query = {'domains': 'derivatives', 'desc': 'preproc',
             'suffix': 'bold', 'extensions': ['.nii', '.nii.gz']}
    if opts.run:
        query['run'] = '|'.join(opts.run)
    if opts.task:
        query['task'] = '|'.join(opts.task)
    if opts.space:
        query['space'] = opts.space
        if opts.space == 'template':
            query['space'] = '|'.join(get_tpl_list())

    # Preprocessed files that are input to the workflow
    prepped_bold = layout.get(**query)
    if not prepped_bold:
        print('No preprocessed files found under the given derivatives '
              'folder "%s".' % derivatives_dir, file=sys.stderr)

    # The magic happens here
    if 'participant' in opts.analysis_level:
        from workflows import first_level_wf
        base_entities = set(['subject', 'session', 'task', 'run', 'acquisition', 'reconstruction'])
        in_files = []
        in_mask = []
        ev_files = []
        regr_files = []
        tr = []

        for part in prepped_bold:
            entities = part.entities
            base = base_entities.intersection(entities)
            subquery = {k: v for k, v in entities.items() if k in base}
            in_files.append(part.path)
            in_mask.append(layout.get(
                domains='derivatives',
                suffix='mask',
                return_type='file',
                space=query['space'],
                **subquery)[0])
            ev_files.append(layout.get(suffix='events', return_type='file', **subquery)[0])
            regr_files.append(layout.get(
                domains='derivatives',
                suffix='regressors',
                return_type='file',
                **subquery)[0])
            tr.append(part.metadata.get('RepetitionTime'))

        workflow = first_level_wf()
        workflow.inputs.inputnode.in_file = in_files[0]
        workflow.inputs.inputnode.in_mask = in_mask[0]
        workflow.inputs.inputnode.regressors = regr_files[0]
        workflow.inputs.inputnode.events_file = ev_files[0]
        workflow.inputs.inputnode.repetition_time = tr[0]
        workflow.base_dir = opts.work_dir
        workflow.run()

    if 'group' in opts.analysis_level:
        pass

    return 0


if __name__ == '__main__':
    sys.exit(main())
