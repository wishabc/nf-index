#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import argparse




def read_masterlist_no_header(masterlist_path):
    dhs_colnames = ['seqname', 'start', 'end', 'id', 'total_signal', 'num_samples', 'num_peaks', 'width', 'summit', 'core_start', 'core_end']
    return pd.read_table(
        masterlist_path,
        header=None, 
        names=dhs_colnames
    )

def main(masterlist, mask, filter_strategy='keep', singleton_cutoff=1):
    if filter_strategy in ('filter_median', 'filter_all'):
        blacklist_filtered = masterlist[mask]
        non_singletons_mask = blacklist_filtered.eval(f'num_samples > {singleton_cutoff}')
        if filter_strategy == 'filter_all':
            valid_dhs = non_singletons_mask
        else:
            cutoff = blacklist_filtered[non_singletons_mask]['mean_signal'].quantile(0.5)
            print(cutoff)
            valid_dhs = blacklist_filtered.eval(f'(num_samples > {singleton_cutoff}) | (mean_signal >= {cutoff})')
        print(np.sum(valid_dhs), valid_dhs.shape)
        mask[mask] = valid_dhs

    masterlist_filtered = masterlist[mask]
    return masterlist_filtered, mask

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter DHS index and create corresponding masks')
    parser.add_argument('masterlist', help='Path to masterlist file')
    parser.add_argument('non_zero_rows_mask', help='Mask marking rows of binary matrix with at least one non-zero element')
    parser.add_argument('blacklist_rows_mask', help='Mask marking rows of binary matrix that are in the blacklist')
    parser.add_argument('filtered_masterlist', help='Resulting filtered masterlist')
    parser.add_argument('filtered_mask', help='Resulting mask to convert masterlist to filtered masterlist')
    parser.add_argument('--singletons_strategy', help='''Strategy for working with singletons in the data.
    One of "keep_all", "filter_median", "filter_all". 
    Filter_median - keep singletons with signal higher than median signal of all reproduced peaks.''', choices=('keep_all', 'filter_median', 'filter_all'), default='keep_all')
    args = parser.parse_args()

    #Read in Rows to Remove for Encode Blacklist
    blacklist_rows = np.loadtxt(args.blacklist_rows_mask, dtype=bool)
    non_zero_rows = np.loadtxt(args.non_zero_rows_mask, dtype=bool)
    valid_rows_mask = ~blacklist_rows & non_zero_rows

    print("Reading Masterlist")
    masterlist = read_masterlist_no_header(args.masterlist)
    masterlist['mean_signal'] = masterlist.eval('total_signal / num_samples')
    print("Filtering Masterlist")
    filtered_index, final_mask = main(masterlist, mask=valid_rows_mask, filter_strategy=args.singletons_strategy)

    filtered_index.to_csv(
        args.filtered_masterlist,
        index=False,
        header=False,
        sep="\t"
    )
    np.savetxt(args.filtered_mask, final_mask, fmt="%d")
