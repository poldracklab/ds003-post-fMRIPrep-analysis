"""
###General functions###
"""

#This function converts fMRIPrep outputs into a Nipype compliant format
#It is in this method that you specify the elements of the fMRIPrep Output that you need for your analysis

# !!!! Note, no modelling is done in this function!
"""
in_file = A pre-processed image file (a single BOLD run)

events_file = a tsv file with event information in BIDS format

confounds_file = a tsv file with regressor information in BIDS derivatives format ([]_desc-confounds_timeseris.tsv)

event_cols = A list of strings that are the column names to be selected from the events_file (do not need to specify onsets, durations or amplitudes, as these are assumed)

regressors_names = A list of strings that are the column names to be selected from the confounds_file. These are *confound* regressors ONLY.

motion_columns = A list of strings that are the column names of motion regressors to be selected from the confounds_file

decimals = An integer specifying rounding

default_amplitude = Default value for event response amplitude
"""

"""relevant events_cols= 'trial_type', 'duration', """
def _bids2nipypeinfo(
    in_file,
    events_file,
    confounds_file,
    events_cols=None,
    regressors_names=None,
    motion_columns=None,
    decimals=3,
    default_amplitude=1.0
):
    from pathlib import Path
    import numpy as np
    import pandas as pd
    from nipype.interfaces.base.support import Bunch

    # Process the events tsv file -- see BIDS
    events = pd.read_csv(events_file, sep='\t', dtype={'trial_type': str})

    events_cols = ['onsets', 'durations', 'amplitudes']
    if events_cols is not None:
        events_cols += events_cols

    # Generates motion column headers for translations and rotations
    if not motion_columns:
        from itertools import product
        motion_columns = ['_'.join(v) for v in product(('trans', 'rot'), 'xyz')]

    # Reads motion.par file and extracts its abso lute path. This is specific for FSL.
    out_motion = Path('motion.par').resolve()

    # TODO: allow for a list of confounds files.
    if isinstance(confounds_file, list):
        confounds_file = confounds_file[0]

    # Process the regressor tsv file -- see BIDS derivatives
    confound_data = pd.read_csv(confounds_file, sep='\t')

    # Pulls motion columns from confounds_file
    np.savetxt(out_motion, confound_data[motion_columns].values, '%g')

    # removing motion columns from the regressor data frame.
    # This may not be necessary, but keeps the data frames clean and avoids redudancy
    # if motion paramaters are treated 'special' by the pipeline this is needed,
    # if they aren't this isn't...
    if regressors_names is None:
        regressors_names = sorted(set(confound_data.columns) - set(motion_columns))

    if regressors_names:
        events_cols += ['regressor_names']
        events_cols += ['regressors']

    # Extract their names and make sure your trial types do not contain spaces
    events['trial_type'] = events.trial_type.fillna('n/a')
    events['trial_type'] = events.trial_type.str.replace(' ', '_')
    condition_names = set(events.trial_type.dropna().values)
    condition_names.discard('n/a')

    # A Bunch is a data structure that allows to access attributes as properties:
    # e.g., runinfo.scans.
    # NiPype requires to be fed Bunches for the modeling interfaces
    runinfo = Bunch(scans=in_file,
                    conditions=sorted(condition_names),
                    **{k: [] for k in events_cols})

    for condition in runinfo.conditions:
        event = events[events.trial_type.str.match(condition)]

        runinfo.onsets.append(np.round(event.onset.values, decimals).tolist())
        runinfo.durations.append(np.round(event.duration.values, decimals).tolist())

        # This will pull in amplitudes *if* they are there.
        if 'amplitudes' in events.columns:
            runinfo.amplitudes.append(np.round(event.amplitudes.values, decimals).tolist())
        else:
            runinfo.amplitudes.append([default_amplitude] * len(event))

    if 'regressor_names' in events_cols:
        runinfo.regressor_names = regressors_names
        try:
            runinfo.regressors = confound_data[regressors_names]
        except KeyError:
            regressors_names = list(set(regressors_names).intersection(
                                    set(confound_data.columns)))
            runinfo.regressors = confound_data[regressors_names]
        runinfo.regressors = runinfo.regressors.fillna(0.0).values.T.tolist()

    return [runinfo], str(out_motion)


def _get_tr(in_dict):
    return in_dict.get('RepetitionTime')


def _len(inlist):
    return len(inlist)


def _dof(inlist):
    return len(inlist) - 1


def _neg(val):
    return -val


def _dict_ds(in_dict, sub, order=['bold', 'mask', 'events', 'regressors', 'tr']):
    return tuple([in_dict[sub][k] for k in order])
