import anndata as ad
import sys
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
import os
import gc


def main(adata, meta, matrices):
    meta['index_peak_files'] = adata.obs['peaks_file']
    adata.obs = meta
    for matrix in matrices:
        matrix_name = os.path.basename(matrix).replace('.npy', '').replace('binary.index.', '')
        matrix = csr_matrix(np.load(matrix).T)
        adata.layers[matrix_name] = matrix

    adata.varm['projected_peaks_binary'] = adata.X.sum(axis=0).A1.squeeze()
    adata.varm['final_qc_passing_dhs'] = (adata.varm['projected_peaks_binary'] > 0) & adata.varm['autosomal_pseudo_reproduced_dhs']
    return adata


if __name__ == '__main__':
    adata_obj = ad.read(sys.argv[1])
    samples_meta = pd.read_table(sys.argv[2]).set_index('ag_id').loc[adata_obj.obs.index]

    matrices = sys.argv[4:]

    result = main(adata_obj, samples_meta, matrices)
    result.write(sys.argv[3])
