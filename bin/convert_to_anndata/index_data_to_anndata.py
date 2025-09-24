import anndata as ad
import sys
import pandas as pd
import numpy as np
from helpers import convert_to_sparse_if_sufficently_sparse, load_from_file
    

def main(rows_meta, cols_meta, matrix):
    adata = ad.AnnData(X=matrix, obs=cols_meta, var=rows_meta)

    adata.obs['n_peaks'] = adata.X.sum(axis=1).A1
    adata.obs['final_qc_passing_sample'] = 1

    adata.var['autosomal_dhs'] = adata.var['#chr'].str.contains('chr[1-9]', regex=True)
    adata.var['n_contributing_samples'] = adata.X.sum(axis=0).A1

    return adata


if __name__ == '__main__':
    annotated_masterlist = pd.read_table(sys.argv[1]).set_index('dhs_id')
    samples_order = np.loadtxt(sys.argv[2], dtype=str)
    
    matrix = load_from_file(sys.argv[3]).T
    matrix = convert_to_sparse_if_sufficently_sparse(matrix)

    peaks_column = sys.argv[5]
    samples_meta = pd.read_table(sys.argv[4]).set_index('ag_id').loc[
        samples_order, [peaks_column]
    ]
    adata_obj = main(annotated_masterlist, samples_meta, matrix)
    adata_obj.write_zarr(sys.argv[6])