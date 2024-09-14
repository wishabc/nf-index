import anndata as ad
import sys
import pandas as pd
import numpy as np
import scipy.sparse as sp
import os


def convert_to_sparse_if_needed(matrix, threshold=0.7):
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
    

def main(rows_meta, cols_meta, matrix, masks):
    adata = ad.AnnData(X=matrix, obs=cols_meta[['peaks_file']], var=rows_meta)
    for mask_name, mask in masks.items():
        if mask_name == 'blacklist_rows':
            mask_name = 'non_blacklisted_dhs'
            mask = ~mask
        elif mask_name == 'filtered.autosomes':
            mask_name = 'autosomal_dhs'
        elif mask_name == 'filtered_DHS':
            mask_name = 'pseudo_reproduced_dhs'
        else:
            continue

        adata.varm[mask_name] = mask
    

    adata.obsm['n_peaks'] = adata.X.sum(axis=1).A1.squeeze()
    adata.varm['n_contributing_samples'] = adata.X.sum(axis=0).A1.squeeze()
    adata.varm['autosomal_pseudo_reproduced_dhs'] = adata.varm['autosomal_dhs'] & adata.varm['pseudo_reproduced_dhs']

    return adata


if __name__ == '__main__':
    annotated_masterlist = pd.read_table(sys.argv[1]).set_index('dhs_id')
    samples_order = np.loadtxt(sys.argv[2], dtype=str)
    matrix = convert_to_sparse_if_needed(np.load(sys.argv[3]).T)

    samples_meta = pd.read_table(sys.argv[4]).set_index('ag_id').loc[samples_order]

    row_masks = {}
    for row_mask in sys.argv[6:]:
        row_mask_name = os.path.basename(row_mask).replace('.mask.txt', '').replace('binary.index.', '')
        mask = np.loadtxt(row_mask, dtype=bool)
        row_masks[row_mask_name] = mask


    adata_obj = main(annotated_masterlist, samples_meta, matrix, row_masks)
    adata_obj.write(sys.argv[5])