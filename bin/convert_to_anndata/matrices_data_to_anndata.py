import sys
import pandas as pd
import os
from helpers import read_zarr_backed, add_matrices_to_anndata


def main(adata, meta, matrices):
    meta['index_peaks_file'] = adata.obs['peaks_file']
    adata.obs = meta
    print(adata)
    matrices_mapping = {
        os.path.basename(matrix).replace('.raw.matrix.npy', ''): matrix 
        for matrix in matrices
        }
    add_matrices_to_anndata(adata, matrices_mapping)


    adata.varm['projected_peaks_binary'] = adata.X.to_memory().sum(axis=0).A1.squeeze()
    adata.varm['final_qc_passing_dhs'] = (adata.varm['projected_peaks_binary'] > 0) & adata.varm['autosomal_pseudo_reproduced_dhs']
    return adata


if __name__ == '__main__':
    adata_obj = read_zarr_backed(sys.argv[1])
    samples_meta = pd.read_table(sys.argv[2]).set_index('ag_id')
    missing_indices = adata_obj.obs.index.difference(samples_meta.index)
    if len(missing_indices) > 0:
        print(f"Warning: Missing indices: {missing_indices}")
    
    samples_meta = samples_meta.loc[adata_obj.obs.index]

    matrices = sys.argv[4:]

    print(adata_obj)

    main(adata_obj, samples_meta, matrices).write_zarr(sys.argv[3])
