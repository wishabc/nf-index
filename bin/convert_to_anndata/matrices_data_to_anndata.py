import sys
import pandas as pd
import os
from helpers import add_matrices_to_anndata
from genome_tools.data.anndata import read_zarr_backed


def main(adata, meta, matrices):
    meta['index_peaks_file'] = adata.obs['peaks_file']
    adata.obs = meta
    matrices_mapping = {
        os.path.basename(matrix).replace('.raw.matrix.npy', ''): matrix 
        for matrix in matrices
        }
    add_matrices_to_anndata(adata, matrices_mapping)


    adata.var['projected_peaks_binary'] = adata.layers['binary'].sum(axis=0).A1.squeeze()
    adata.var['final_qc_passing_dhs'] = (adata.var['projected_peaks_binary'] > 0) & adata.var['autosomal_dhs']
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
