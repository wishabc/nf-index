
import zarr
import anndata as ad
import dask.array as da
import numpy as np
import scipy.sparse as sp


def _read_zarr_group(dat):
    if type(dat) == zarr.Array:
        return da.from_zarr(dat)
    else:
        return ad.experimental.sparse_dataset(dat)


def read_zarr_backed(filename):
    d = {}
    attributes = ["obsm", "varm", "obsp", "varp", "uns"]
    df_attributes = ["obs", "var"]

    with zarr.open(filename, "r") as f:
        if "encoding-type" in f.attrs:
            attributes.extend(df_attributes)
        
        d['X'] = _read_zarr_group(f['X'])
        for k in attributes:
            if k not in f:
                continue
            d[k] = ad.experimental.read_elem(f[k])

        if 'layers' in f.keys():
            d['layers'] = {}
            for layer in f['layers']:
                d['layers'][layer] = _read_zarr_group(f['layers'][layer])

        adata = ad.AnnData(**d)

    return adata


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
        if matrix.shape[0] != adata.shape[0]:
            raise ValueError(f"Matrix {matrix_name} has {matrix.shape[0]} samples, but the number of samples in anndata is {adata.shape[0]}")
        if matrix.shape[1] != adata.shape[1]:
            mask = adata.varm['final_qc_passing_dhs']
            assert mask.sum() == matrix.shape[1]
            data = np.empty(adata.shape, dtype=matrix.dtype)
            data[:] = np.nan
            data[:, mask] = matrix
            matrix = data

        adata.layers[matrix_name] = convert_to_sparse_if_sufficently_sparse(matrix)
