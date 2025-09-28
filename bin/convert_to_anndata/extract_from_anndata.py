import scipy.sparse as sp
import numpy as np
import argparse
import pandas as pd
import dask.array as da
from genome_tools.data.anndata import read_zarr_backed
from helpers import get_mask_from_column_name

def main(anndata_obj, extra_keys):
    metadata = anndata_obj.obs
    masterlist = anndata_obj.var.reset_index()[['#chr', 'start', 'end', 'dhs_id', 'dhs_summit']]
    
    matrices = {}
    for key in extra_keys:
        if key in anndata_obj.layers:
            matrices[key] = anndata_obj.layers[key]
        elif key == "bg_corrected_agg_cutcounts":
            val = anndata_obj.obs['nuclear_reads'].values[:, None] * anndata_obj.layers['density'].compute() / 1e6 - anndata_obj.layers['mean_bg_agg_cutcounts'].compute()
            val = np.nan_to_num(val, copy=False)
            np.clip(val, 0, None, out=val)
            matrices[key] = val
        else:
            raise ValueError(f"Layer '{key}' not found in the AnnData object.")

    return metadata, masterlist, matrices


def convert_matrix_to_dense(matrix):
    if isinstance(matrix, da.Array):
        # If the matrix is a Dask array, compute it to convert to a NumPy array
        return matrix.compute()
    elif sp.issparse(matrix):
        # If the matrix is a Scipy sparse matrix, convert to a dense NumPy array
        return matrix.toarray()
    elif isinstance(matrix, np.ndarray):
        # If the matrix is already a NumPy array, return it as is
        return matrix
    else:
        raise TypeError("Unsupported matrix type. Must be a Dask array, NumPy array, or Scipy sparse matrix.")


if __name__ == '__main__':
    args = argparse.ArgumentParser('Extract data from anndata object containing DHS information')
    args.add_argument('anndata', help='Path to anndata zarr object')
    args.add_argument('index', help='Path where to save DHS index')
    args.add_argument('samples_order', help='Path where to save samples order')
    args.add_argument('samples_meta', help='Path where to save samples metadata')
    args.add_argument('--extra_layers', nargs='+', help='Names of extra layers to extract from anndata', default=[])
    args.add_argument('--extra_layers_suffix', help='Suffix to add to the extra layers names', default='matrix')
    args.add_argument('--dhs_mask_name', help='Name of the var layer containing DHS mask', default="")
    args.add_argument('--matrix_samples_file', help='Name of the var layer containing DHS mask', default=None)
    args = args.parse_args()
    anndata = read_zarr_backed(args.anndata)
    mask = get_mask_from_column_name(anndata, args.dhs_mask_name)

    anndata = anndata[:, mask]
    
    metadata, index, matrices = main(anndata, args.extra_layers)
    index.to_csv(args.index, sep='\t', index=False, header=False)
    if args.matrix_samples_file is not None:
        matrix_metadata = pd.read_table(args.matrix_samples_file).set_index('ag_id')
        missing_ids = matrix_metadata.index.difference(metadata.index)
        print(f"Missing IDs in matrix samples metadata: {missing_ids}")
        if matrix_metadata.shape[0] == metadata.shape[0] and len(missing_ids) == 0:
            pass
        else:
            if len(missing_ids) > 0:
                print(f"Missing IDs in anndata samples metadata: {missing_ids}")
            else:
                print("Samples metadata does not match the number of samples in the anndata object")
            print("Ignoring anndata samples metadata")
            metadata = matrix_metadata

    np.savetxt(args.samples_order, metadata.index, fmt='%s')
    metadata.reset_index().to_csv(args.samples_meta, sep='\t', index=False)

    for name, matrix in matrices.items():
        if sp.issparse(matrix):
            sp.save_npz(f"{name}.{args.extra_layers_suffix}.npz", matrix.T)
        elif isinstance(matrix, da.Array):
            np.save(f"{name}.{args.extra_layers_suffix}.npy", matrix.T.compute())
        else:
            np.save(f"{name}.{args.extra_layers_suffix}.npy", matrix.T)
        