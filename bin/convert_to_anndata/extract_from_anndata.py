import sys
import anndata as ad
import numpy as np
import argparse

def main(anndata_obj, extra_keys):
    samples_order = anndata_obj.obs.index
    masterlist = anndata_obj.var.reset_index()[['#chr', 'start', 'end', 'dhs_id']]
    
    matrices = {key: anndata_obj.layers[key] for key in extra_keys}
    return samples_order, masterlist, matrices


if __name__ == '__main__':
    args = argparse.ArgumentParser('Extract data from anndata object containing DHS information')
    args.add_argument('anndata', help='Path to anndata object')
    args.add_argument('index', help='Path where to save DHS index')
    args.add_argument('samples_order', help='Path where to save samples order')
    args.add_argument('--extra_layers', nargs='+', help='Names of extra layers to extract from anndata', default=[])
    args.add_argument('--extra_layers_suffix', help='Suffix to add to the extra layers names', default='matrix')
    args.add_argument('--dhs_mask_name', help='Name of the varm layer containing DHS mask')
    args = args.parse_args()
    anndata = ad.read(args.anndata)
    anndata = anndata[:, anndata.varm[args.dhs_mask_name]]
    
    samples_order, index, matrices = main(anndata, args.extra_layers)
    index.to_csv(args.index, sep='\t', index=False, header=False)
    np.savetxt(args.samples_order, samples_order, fmt='%s')

    for name, matrix in matrices.items():
        np.save(f"{name}.{args.extra_layers_suffix}.npy", matrix.T)