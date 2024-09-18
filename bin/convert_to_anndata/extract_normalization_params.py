import numpy as np
import sys
from helpers import read_zarr_backed
import base64
import json


def main(anndata_obj):
    deseq_params = base64.b64decode(anndata_obj.uns['deseq_params'])
    
    norm_json = anndata_obj.uns['lowess_params']
    norm_arrays = {}
    for key in anndata_obj.uns.keys():
        if key.startswith('norm_params_'):
            real_key = key.replace('norm_params_', '')
            norm_arrays[real_key] = anndata_obj.uns[key]

    return deseq_params, norm_arrays, norm_json



if __name__ == '__main__':
    
    anndata = read_zarr_backed(sys.argv[1])
    savedir = sys.argv[2]
    
    deseq_params, norm_arrays, norm_json = main(anndata)
    with open(f'{savedir}/deseq_params.RDS', 'wb') as f:
        f.write(deseq_params)
    
    with open(f'{savedir}/lowess_params.json', 'w') as f:
        json.dump(norm_json, f)
    
    np.savez_compressed(f'{savedir}/lowess_params.npz', **norm_arrays)