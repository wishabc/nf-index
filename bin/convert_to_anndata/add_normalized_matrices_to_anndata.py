import sys
import pandas as pd
import numpy as np
import os
from helpers import read_zarr_backed, add_matrices_to_anndata
import json


def main(adata, matrices, params, formula, annotated_masterlist):
    matrices_mapping = {
        os.path.basename(matrix).replace(
            'normalized.only_autosomes.filtered.', ''
        ).replace(
            '.npy', ''
        ): matrix for matrix in matrices
    }
    add_matrices_to_anndata(adata, matrices_mapping)
    for param in params:
        if param.endswith('RDS'):
            with open(param, 'rb') as f:
                adata.uns['deseq_params'] = f.read()
        elif param.endswith('json'):
            with open(param) as f:
                adata.uns['lowess_params'] = json.load(f)
        elif param.endswith('npz'):
            loaded_params = np.load(param)
            for key in loaded_params:
                adata.uns[f'norm_params_{key}'] = loaded_params[key]
        else:
            raise ValueError(f'Unknown parameter file type: {param}')
    adata.uns['vp_formula'] = formula
    adata.var = adata.var.join(annotated_masterlist, how='left')
    return adata


if __name__ == '__main__':
    adata = read_zarr_backed(sys.argv[1])
    annotated_masterlist = pd.read_table(sys.argv[3]).set_index('V4')
    annotated_masterlist = annotated_masterlist[
        [x for x in annotated_masterlist.columns if not x.startswith('V')]
    ]
    annotated_masterlist.rename(columns={x: f'{x}_variance_partition' for x in annotated_masterlist.columns}, inplace=True)
    formula = sys.argv[4]
    matrices = sys.argv[5:8]
    params = sys.argv[8:10]

    adata = main(adata, matrices, params, formula, annotated_masterlist)
    adata.write_zarr(sys.argv[2])