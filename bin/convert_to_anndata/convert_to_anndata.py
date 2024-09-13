import anndata as ad
import sys
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
import os


def main(rows_meta, cols_meta, matrix, masks):
    adata = ad.AnnData(X=matrix.T, obs=cols_meta, var=rows_meta)
    for mask_name, mask in masks.items():
        adata.varm[mask_name] = mask

    return adata


if __name__ == '__main__':
    annotated_masterlist = pd.read_table(sys.argv[1]).set_index('dhs_id')
    samples_order = np.loadtxt(sys.argv[2], dtype=str)
    #matrix = csr_matrix(np.load(sys.argv[3]))
    matrix = np.zeros((len(annotated_masterlist), len(samples_order)), dtype=bool)

    samples_meta = pd.read_table(sys.argv[4]).set_index('ag_id').loc[samples_order]

    outpath = sys.argv[5]

    row_masks = {}
    for row_mask in sys.argv[6:]:
        row_mask_name = os.path.basename(row_mask).replace('.mask.txt', '')
        mask = np.loadtxt(row_mask, dtype=bool)
        row_masks[row_mask_name] = mask


    adata_obj = main(annotated_masterlist, samples_meta, matrix, row_masks)
    adata_obj.write(outpath)