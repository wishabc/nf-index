#!/usr/bin/env python3

#Load Dependencies
import pandas as pd
import numpy as np
import sys

masterlist_prefix = sys.argv[1]
percentile = sys.argv[2]

#Read in Rows to Remove for Encode Blacklist
blacklist_rows = pd.read_table(sys.argv[3], header=None, names=["row_id"])
non_zero_rows = pd.read_table(sys.argv[4], header=None, names=["row_id"])
blacklist_rows = blacklist_rows * non_zero_rows
#Read in masterlist file for calculating average score
print("Reading Masterlist")
masterlist = pd.read_table(sys.argv[5], header=None, names=['seqname', 'start', 'end', 'id', 'total_signal', 'num_samples', 'num_peaks', 'width', 'summit', 'core_start', 'core_end'])
print(masterlist.shape[0])
m_size = masterlist.shape[0]

#Remove blacklist rows from dataframe but keep index
masterlist['index_column'] = masterlist.index
blacklist_filtered = masterlist.loc[~masterlist['index_column'].isin(blacklist_rows['row_id'])]
blacklist_filtered = blacklist_filtered[['seqname', 'start', 'end', 'id', 'total_signal', 'num_samples', 'num_peaks', 'width', 'summit', 'core_start', 'core_end']]
blacklist_filtered['mean_signal'] = blacklist_filtered['total_signal'] / blacklist_filtered['num_samples']
blacklist_filtered.to_csv(str(masterlist_prefix) + '.blacklistfiltered.bed', index=False, header=False, sep="\t")
blacklist_filtered = blacklist_filtered[['total_signal', 'num_samples', 'mean_signal']]

#Calculate cutoff score for threshold
#Cutoff score = score where ${ercentile} of peaks replicate
i = blacklist_filtered['num_samples'] >= 2
cutoff = blacklist_filtered[i]['mean_signal'].quantile(float(percentile))
print(percentile)
print(cutoff)

#How many DHS's are greater than the threshold
j= (blacklist_filtered['num_samples'] >= 2) | ((blacklist_filtered['num_samples']== 1) & (blacklist_filtered['mean_signal']>=cutoff))
print(np.sum(j))
print(j)

#Print Filtered Masterlist
index_subset = blacklist_filtered[j]
k = index_subset.index
masterlist_filtered = masterlist.loc[k]
masterlist_filtered['mean_signal'] = masterlist_filtered['total_signal'] / masterlist_filtered['num_samples']
new_masterlist = masterlist_filtered[['seqname', 'start', 'end', 'id', 'total_signal', 'num_samples', 'num_peaks', 'width', 'summit', 'core_start', 'core_end', 'mean_signal']]
new_masterlist.to_csv(sys.argv[6], index=False, header=False, sep="\t")

#filtered_binary = binary.loc[binary['index_column'].isin(new_masterlist.index)]
#filtered_binary = filtered_binary.reset_index(drop=True)
#filtered_binary = filtered_binary.drop(columns=['index_column'])

s = pd.DataFrame(k, columns=['row'])
d = pd.DataFrame(np.zeros((m_size, 1)), columns=['all'])
d.loc[s['row'], 'all'] = 1
d['all'] = d['all'].astype(int)
d.to_csv(sys.argv[7], index=False, header=False, sep="\t")
