import numpy as np
import scipy.sparse as sp


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
    

def add_matrices_to_anndata(adata, matrices_mapping):
    for matrix_name, matrix in matrices_mapping.items():
        matrix = np.load(matrix).T
        if matrix.dtype == np.float64:
            matrix = matrix.astype(np.float32)

        if matrix.shape[0] != adata.shape[0]:
            raise ValueError(f"Matrix {matrix_name} has {matrix.shape[0]} samples, but the number of samples in anndata is {adata.shape[0]}")
        if matrix.shape[1] != adata.shape[1]:
            mask = adata.var['final_qc_passing_dhs']
            assert mask.sum() == matrix.shape[1]
            data = np.empty(adata.shape, dtype=np.float32)
            data[:] = np.nan
            data[:, mask] = matrix
            matrix = data

        adata.layers[matrix_name] = convert_to_sparse_if_sufficently_sparse(matrix)
