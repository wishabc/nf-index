import numpy as np
import sys
import base64
import json
from genome_tools.data.anndata import read_zarr_backed


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.bool_):  # Handling numpy.bool_
            return bool(obj)
        return super(NpEncoder, self).default(obj)
    

def main(anndata_obj):
    deseq_params = base64.b64decode(anndata_obj.uns['vst_dispersion_function_rds.b64_encoded'])
    
    norm_json = anndata_obj.uns['lowess_normalization_params']
    norm_arrays = {}
    for key in anndata_obj.uns.keys():
        if key.startswith('lowess_normalization_params.'):
            real_key = key.replace('lowess_normalization_params.', '')
            norm_arrays[real_key] = anndata_obj.uns[key]

    return deseq_params, norm_arrays, norm_json



if __name__ == '__main__':
    
    anndata = read_zarr_backed(sys.argv[1])
    savedir = sys.argv[2]
    
    deseq_params, norm_arrays, norm_json = main(anndata)
    with open(f'{savedir}/deseq_params.RDS', 'wb') as f:
        f.write(deseq_params)
    
    with open(f'{savedir}/lowess_params.json', 'w') as f:
        json.dump(norm_json, f, cls=NpEncoder, indent=2)
    
    np.savez_compressed(f'{savedir}/lowess_params.npz', **norm_arrays)