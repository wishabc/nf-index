import sys
import anndata as ad
import numpy as np


def main(anndata_obj, extra_keys):
    samples_order = anndata_obj.obs.index
    masterlist = anndata_obj.var.reset_index()[['#chr', 'start', 'end', 'dhs_id']]
    
    matrices = {key: anndata_obj.layers[key] for key in extra_keys}
    return samples_order, masterlist, matrices

# FIXME add argparse
if __name__ == '__main__':
    anndata = ad.read(sys.argv[1])
    anndata = anndata[:, anndata.varm[sys.argv[4]]]
    
    extra_keys = sys.argv[5:]
    samples_order, index, matrices = main(anndata, extra_keys)
    index.to_csv(sys.argv[2], sep='\t', index=False, header=False)
    np.savetxt(sys.argv[3], samples_order, fmt='%s')
    for name, matrix in matrices.items():
        np.save(f"{name}.matrix.npy", matrix.T)