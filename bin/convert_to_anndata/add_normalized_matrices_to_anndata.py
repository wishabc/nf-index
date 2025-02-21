import sys
import pandas as pd
import numpy as np
import os
from helpers import add_matrices_to_anndata
from genome_tools.data.anndata import read_zarr_backed
import json
import base64

def main(adata, matrices, params, formula, annotated_masterlist, mask_name='autosomes_mask'):
    matrices_mapping = {
        os.path.basename(matrix).replace(
            'normalized.only_autosomes.filtered.', ''
        ).replace(
            '.npy', ''
        ): matrix for matrix in matrices
    }
    add_matrices_to_anndata(adata, matrices_mapping, mask_name)
    for param in params:
        if param.endswith('RDS'):
            with open(param, 'rb') as f:
                adata.uns['deseq_params'] = base64.b64encode(f.read()).decode('utf-8')
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
    assert len(annotated_masterlist) == len(adata.var)
    annotated_masterlist = annotated_masterlist[
        [x for x in annotated_masterlist.columns if not x.startswith('V')]
    ]
    annotated_masterlist.rename(columns={x: f'{x}_variance_partition' for x in annotated_masterlist.columns}, inplace=True)
    formula = sys.argv[4]
    mask_name = sys.argv[5]
    matrices = sys.argv[6:9]
    params = sys.argv[9:]

    adata = main(adata, matrices, params, formula, annotated_masterlist, mask_name=mask_name)
    adata.write_zarr(sys.argv[2])