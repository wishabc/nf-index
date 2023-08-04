#!/usr/bin/env python3

#Load Dependencies
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

masterlist_prefix = sys.argv[1]
percentile = sys.argv[2]

#Read in Rows to Remove for Encode Blacklist
blacklist_rows = pd.read_table('blacklist_rows.txt', header=None)

#Read in masterlist file for calculating average score
print("Reading Masterlist")
masterlist = pd.read_table('masterlist_DHSs_' + str(masterlist_prefix) + '_all_chunkIDs.bed', header=None)
print(masterlist.shape[0])

#Assign headers to dataframes
blacklist_rows.columns= ["row_id"]
masterlist.columns = ['seqname', 'start', 'end', 'id', 'total_signal', 'num_samples', 'num_peaks', 'width', 'summit', 'core_start', 'core_end']

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

#Print Filtered Masterlist
index_subset = blacklist_filtered[j]
k = index_subset.index
print(k)
masterlist_filtered = masterlist.loc[k]
masterlist_filtered['mean_signal'] = masterlist_filtered['total_signal'] / masterlist_filtered['num_samples']
new_masterlist = masterlist_filtered[['seqname', 'start', 'end', 'id', 'total_signal', 'num_samples', 'num_peaks', 'width', 'summit', 'core_start', 'core_end', 'mean_signal']]

#Print Filtered Binary Matrix
#filtered_binary = binary.loc[binary['index_column'].isin(new_masterlist.index)]
#filtered_binary = filtered_binary.reset_index(drop=True)
#filtered_binary = filtered_binary.drop(columns=['index_column'])

new_masterlist.to_csv('masterlist_DHSs_' + str(masterlist_prefix) + str(percentile) + '.blacklistfiltered.bed', index=False, header=False, sep="\t")
#filtered_binary.to_csv(str(masterlist_prefix) + str(percentile) + '.blacklistfiltered.binary.mtx', index=False, header=False, sep="\t")

def create_mask(arr):
    result = []
    for i in range(len(arr) - 1):
        result.append(arr[i])
        if arr[i] + 1 != arr[i + 1]:
            # Fill the gap with zeros
            gap_size = arr[i + 1] - arr[i] - 1
            result.extend([0] * gap_size)
	else:
	   result.append(1)
    result.append(arr[-1])  # Add the last element of the original array
    return result

mask = create_mask(k)
np.savetxt(fname="masked_elements", mask, delimiter='\t')
