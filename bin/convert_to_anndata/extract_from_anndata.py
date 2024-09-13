import sys
import anndata as ad
import numpy as np


def main(anndata_obj):
    samples_order = anndata_obj.obs.index
    masterlist = anndata_obj.var.reset_index()[['#chr', 'start', 'end', 'dhs_id']]
    return samples_order, masterlist


if __name__ == '__main__':
    anndata = ad.read(sys.argv[1])
    samples_order, index = main(anndata)
    index.to_csv(sys.argv[2], sep='\t', index=False, header=False)
    np.savetxt(sys.argv[3], samples_order, fmt='%s')