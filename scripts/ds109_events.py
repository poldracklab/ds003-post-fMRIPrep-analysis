from pathlib import Path
from numpy import nan
import pandas as pd

if __name__ == '__main__':
    from argparse import ArgumentParser
    from argparse import RawTextHelpFormatter

    parser = ArgumentParser(description='DS000109 Analysis Workflow',
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument('in_file', action='store', type=Path,
                        help='the input events file we want to fix')

    opts = parser.parse_args()

    original_df = pd.read_csv(opts.in_file, sep='\t', na_values="n/a",
                              dtype={'onset': float, 'duration': float,
                                     'trial_type': str})

    opts.in_file.rename('%s.bak' % opts.in_file)

    tt_names = original_df[~original_df.trial_type.isnull()]
    working_df = original_df[original_df.trial_type.isnull()].set_index('onset')

    for tt in set(tt_names.trial_type.values):
        rows = tt_names[tt_names.trial_type == tt]
        for _, row in rows.iterrows():
                onset = row.onset
                end = (onset + row.duration) - 2.0
                working_df.loc[onset:end, 'trial_type'] = tt

    working_df[working_df.trial_type == 0.0] = nan

    working_df.to_csv(opts.in_file, sep='\t', na_rep='n/a')
