import sys
import pandas as pd
import os
from helpers import add_matrices_to_anndata
from genome_tools.data.anndata import read_zarr_backed
import scipy.sparse as sp
import anndata as ad


def main(adata, matrices):
    matrices_mapping = {
        os.path.basename(matrix).replace('matrix.', '').replace('.npy', ''): matrix 
        for matrix in matrices
    }
    add_matrices_to_anndata(adata, matrices_mapping)

    sums = adata.layers['binary'].sum(axis=0)
    if sp.issparse(adata.layers['binary']):
        sums = sums.A
    adata.var['projected_peaks_binary'] = sums.squeeze()
    adata.var['final_qc_passing_dhs'] = (adata.var['projected_peaks_binary'] > 0) & adata.var['autosomal_dhs']
    return adata


if __name__ == '__main__':
    adata_obj = read_zarr_backed(sys.argv[1])
    samples_meta = pd.read_table(sys.argv[2], low_memory=False).set_index('ag_id')
    samples_meta[
        samples_meta.select_dtypes(include=['object', 'O', 'category']).columns
    ] = samples_meta.select_dtypes(include=['object', 'O', 'category']).fillna('None')
    
    
    missing_indices = adata_obj.obs.index.difference(samples_meta.index)
    if len(missing_indices) > 0:
        print(f"Warning: Missing indices: {missing_indices}")
        print("Creating a new anndata object (dropping X and obs)")
        adata_obj = ad.AnnData(X=None, obs=samples_meta, var=adata_obj.var)
    else:
        samples_meta['peaks_for_index'] = adata_obj.obs['peaks_for_index']
        samples_meta['n_peaks'] = adata_obj.obs['n_peaks']
        adata_obj.obs = samples_meta

    matrices = sys.argv[4:]

    main(adata_obj, matrices).write_zarr(sys.argv[3])
