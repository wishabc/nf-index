#!/usr/bin/env python3

#Load Dependencies
import pandas as pd
import numpy as np
import sys

masterlist_prefix = sys.argv[1]
percentile = float(sys.argv[2])

#Read in Rows to Remove for Encode Blacklist
blacklist_rows = np.loadtxt(sys.argv[3], dtype=bool)
non_zero_rows = np.loadtxt(sys.argv[4], dtype=bool)
valid_rows_mask = ~blacklist_rows & non_zero_rows

print("Reading Masterlist")
masterlist = pd.read_table(sys.argv[5], header=None, names=['seqname', 'start', 'end', 'id', 'total_signal', 'num_samples', 'num_peaks', 'width', 'summit', 'core_start', 'core_end'])
print(masterlist.shape[0])
blacklist_filtered = masterlist[valid_rows_mask]

blacklist_filtered = blacklist_filtered[['seqname', 'start', 'end', 'id', 'total_signal', 'num_samples', 'num_peaks', 'width', 'summit', 'core_start', 'core_end']]
blacklist_filtered['mean_signal'] = blacklist_filtered['total_signal'] / blacklist_filtered['num_samples']
blacklist_filtered.to_csv(f'{masterlist_prefix}.blacklistfiltered.bed', index=False, header=False, sep="\t")
blacklist_filtered = blacklist_filtered[['total_signal', 'num_samples', 'mean_signal']]

singleton_cutoff = 1
#Calculate cutoff score for threshold
#Cutoff score = score where ${ercentile} of peaks replicate
i = blacklist_filtered.eval(f'num_samples > {singleton_cutoff}')
cutoff = blacklist_filtered[i]['mean_signal'].quantile(percentile)
print(percentile)
print(cutoff)

#How many DHS's are greater than the threshold
valid_dhs = blacklist_filtered.eval(f'(num_samples > {singleton_cutoff}) | (num_samples <= {singleton_cutoff} & mean_signal >= {cutoff})')
print(np.sum(valid_dhs))
print(valid_dhs.shape)

valid_rows_mask[valid_rows_mask] = valid_dhs
#Print Filtered Masterlist
masterlist_filtered = masterlist[valid_rows_mask]
masterlist_filtered['mean_signal'] = masterlist_filtered.eval('total_signal / num_samples')
masterlist_filtered[['seqname', 'start', 'end', 'id', 'total_signal', 'num_samples', 'num_peaks', 'width', 'summit', 'core_start', 'core_end', 'mean_signal']].to_csv(sys.argv[6], index=False, header=False, sep="\t")

np.savetxt(sys.argv[7], valid_rows_mask, fmt="%d")
