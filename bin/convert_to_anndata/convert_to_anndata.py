import anndata as ad
import sys
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
import os


def main(rows_meta, cols_meta, matrix, masks):
    adata = ad.AnnData(X=matrix, obs=cols_meta[['peaks_file']], var=rows_meta)
    for mask_name, mask in masks.items():
        if mask_name == 'blacklist_rows':
            mask_name = 'non_blacklisted_dhs'
            mask = ~mask
        elif mask_name == 'filtered.autosomes':
            mask_name = 'autosomal_dhs'
        elif mask_name == 'filtered_DHS':
            mask_name = 'pseudo_reproduced_dhs'
        else:
            continue

        adata.varm[mask_name] = mask
    
    adata.obsm['n_peaks'] = adata.X.tocsc().sum(axis=0)
    adata.varm['n_contributing_samples'] = adata.X.sum(axis=0)
    adata.varm['autosomal_pseudo_reproduced_dhs'] = adata.varm['autosomal_dhs'] & adata.varm['pseudo_reproduced_dhs']

    return adata


if __name__ == '__main__':
    annotated_masterlist = pd.read_table(sys.argv[1]).set_index('dhs_id')
    samples_order = np.loadtxt(sys.argv[2], dtype=str)
    matrix = csr_matrix(np.load(sys.argv[3]).T)

    samples_meta = pd.read_table(sys.argv[4]).set_index('ag_id').loc[samples_order]

    row_masks = {}
    for row_mask in sys.argv[6:]:
        row_mask_name = os.path.basename(row_mask).replace('.mask.txt', '').replace('binary.index.', '')
        mask = np.loadtxt(row_mask, dtype=bool)
        row_masks[row_mask_name] = mask


    adata_obj = main(annotated_masterlist, samples_meta, matrix, row_masks)
    adata_obj.write(sys.argv[5])