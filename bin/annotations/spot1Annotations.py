#!/usr/bin/env python3

#Load Dependencies
import pandas as pd
import numpy as np
import gzip
import sys
import gc
#Load in Binary Matrix and Metadata for subsampled SPOT score

b = np.load(sys.argv[1]).astype(float)

mask_file = sys.argv[2]
mask_df = pd.read_table(mask_file, header=None, names=['mask'])
#b_filtered = b[mask_df['mask'] == 1]
print(b.shape)

samples_file = sys.argv[3]
metadata = pd.read_table(samples_file)

#Make sure the spot1 scores and matrix are matrix compatible
spot1 = metadata['SPOT1_score'].to_numpy()
b *= spot1[None, :]

#Replace 0's with nan's for computation
b = np.where(b != 0, b, np.nan)

#Calculate metrics across rows
spot1_std = np.nanstd(b, axis=1)
spot1_mean = np.nanmean(b, axis=1)

spot1_q = np.nanpercentile(b, [0, 0.25, 0.5, 0.75, 1], axis=1)

del b
gc.collect()
# Combine all metrics in a dataframe

order = ['spot1_std', 'spot1_min', 'spot1_mean', 'spot1_median', 'spot1_max', 'spot1_Q1', 'spot1_Q3']

spot1_metrics = pd.DataFrame({
    'spot1_std': spot1_std,
    'spot1_min': spot1_q[0],
    'spot1_mean': spot1_mean,
    'spot1_median': spot1_q[2],
    'spot1_max': spot1_q[4],
    'spot1_Q1': spot1_q[1],
    'spot1_Q3': spot1_q[3]
})

#Output Results
spot1_metrics[order].to_csv("spot1_metrics.tsv", index=False, header=False, sep="\t")
