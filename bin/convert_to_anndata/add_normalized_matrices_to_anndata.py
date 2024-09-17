import sys
import pandas as pd
import numpy as np
import os
from helpers import read_zarr_backed, add_matrices_to_anndata
import rpy2.robjects as ro
from rpy2.robjects import r
import json


def main(adata, matrices, params):

    adata.shape = adata.varm['final_qc_passing_dhs']
    matrices_mapping = {
        os.path.basename(matrix).replace(
            '.only_autosomes.filtered', ''
        ).replace(
            '.npy', ''
        ): matrix for matrix in matrices
    }
    
    add_matrices_to_anndata(adata, matrices_mapping)
    for param in params:
        if param.endswith('RDS'):
            #adata.uns['deseq_params'] = r.readRDS(param)
            adata.uns['deseq_params'] = param
        elif param.endswith('json'):
            with open(param) as f:
                adata.uns['lowess_params'] = json.load(f)
        elif param.endswith('npy'):
            adata.uns['norm_values'] = np.load(param)
        else:
            raise ValueError(f'Unknown parameter file type: {param}')

    return adata


if __name__ == '__main__':
    adata = read_zarr_backed(sys.argv[1])
    matrices = sys.argv[2:5]
    params = sys.argv[5:7]
    adata = main(adata, matrices, params)
    adata.write_zarr(sys.argv[8])