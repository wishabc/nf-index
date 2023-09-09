#!/usr/bin/env python3

#Load Dependencies
import pandas as pd
import numpy as np
import gzip
import sys

#Load in Binary Matrix and Metadata for subsampled SPOT score
binary_file = sys.argv[1]
with gzip.open(binary_file, 'rb') as f:
    b = np.genfromtxt(f, delimiter='\t', dtype=np.uint8)


samples_file = sys.argv[2]
metadata = pd.read_table(samples_file)

#Make sure the spot1 scores and matrix are matrix compatible
spot1 = metadata['subsampled_hotspot1'].values.reshape(1, -1)
values = b * spot1

#Replace 0's with nan's for computation
values_nonzero = np.where(values != 0, values, np.nan)

#Calculate metrics across rows
spot1_std = np.nanstd(values_nonzero, axis=1)
spot1_min = np.nanmin(values_nonzero, axis=1)
spot1_mean = np.nanmean(values_nonzero, axis=1)
spot1_median = np.nanmedian(values_nonzero, axis=1)
spot1_max = np.nanmax(values_nonzero, axis=1)
spot1_25 = np.nanpercentile(values_nonzero, 25, axis=1)
spot1_75 = np.nanpercentile(values_nonzero, 75, axis=1)

# Combine all metrics in a dataframe
spot1_df = {
    'spot1_std': spot1_std,
    'spot1_min': spot1_min,
    'spot1_mean': spot1_mean,
    'spot1_median': spot1_median,
    'spot1_max': spot1_max,
    'spot1_Q1': spot1_25,
    'spot1_Q3': spot1_75
}

spot1_metrics = pd.DataFrame(spot1_df)

#Output Results
spot1_metrics.to_csv("spot1_metrics.tsv", index=False, header=False, sep="\t")
