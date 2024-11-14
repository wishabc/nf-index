#!/usr/bin/env python3

#Load Dependencies
import pandas as pd
import numpy as np
import sys
import scipy.sparse as sp
import gc


# TODO: move to genometools
def convert_to_sparse_if_sufficently_sparse(matrix, threshold=0.7):
    """
    Check the sparsity of a matrix and convert to sparse format if the proportion of zero values exceeds the threshold.

    Parameters:
    -----------
    matrix: np.ndarray or scipy.sparse.spmatrix
        Input matrix to check for sparsity.
    threshold: float, optional (default=0.y)
        The threshold for sparsity. If the proportion of zero values exceeds this threshold, 
        the matrix will be converted to sparse format.

    Returns:
    --------
    matrix: np.ndarray or scipy.sparse.spmatrix
        The original matrix or the converted sparse matrix.
    """
    # Check if the matrix is already sparse
    if sp.issparse(matrix):
        print("Matrix is already sparse.")
        return matrix
    
    # Calculate the proportion of zeros
    zero_count = np.count_nonzero(matrix == 0)
    total_elements = matrix.size
    sparsity = zero_count / total_elements

    print(f"Sparsity: {sparsity:.2f}")

    # Convert to sparse format if sparsity exceeds the threshold
    if sparsity > threshold:
        print("Converting to sparse format.")
        return sp.csr_matrix(matrix)
    else:
        print("Matrix is sufficiently dense; keeping the current format.")
        return matrix


#Load in Binary Matrix and Metadata for subsampled SPOT score
b = np.load(sys.argv[1])
b = convert_to_sparse_if_sufficently_sparse(b).astype(float)
samples_order = np.loadtxt(sys.argv[2], dtype=str)

metadata = pd.read_table(sys.argv[3]).set_index('ag_id').loc[samples_order]

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
spot1_metrics[order].to_csv(sys.argv[4], index=False, sep="\t")
