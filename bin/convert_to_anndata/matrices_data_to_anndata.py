import sys
import pandas as pd
import numpy as np
import os
from helpers import read_zarr_backed, convert_to_sparse_if_sufficently_sparse


def main(adata, meta, matrices):
    meta['index_peak_files'] = adata.obs['peaks_file']
    adata.obs = meta
    for matrix in matrices:
        matrix_name = os.path.basename(matrix).replace('.raw.matrix.npy', '')
        matrix = convert_to_sparse_if_sufficently_sparse(np.load(matrix).T)
        adata.layers[matrix_name] = matrix

    adata.varm['projected_peaks_binary'] = adata.X.to_memory().sum(axis=0).A1.squeeze()
    adata.varm['final_qc_passing_dhs'] = (adata.varm['projected_peaks_binary'] > 0) & adata.varm['autosomal_pseudo_reproduced_dhs']
    return adata


if __name__ == '__main__':
    adata_obj = read_zarr_backed(sys.argv[1])
    samples_meta = pd.read_table(sys.argv[2]).set_index('ag_id').loc[adata_obj.obs.index]

    matrices = sys.argv[4:]

    main(adata_obj, samples_meta, matrices).write_zarr(sys.argv[3])
