import sys
import numpy as np
import pandas as pd
from genome_tools.data.anndata import read_zarr_backed
from statsmodels.stats.multitest import multipletests


def acat_equal(p, ones_mask):
    denom = (~ones_mask).sum(axis=0)
    num = np.sum(np.tan((0.5 - p) * np.pi), axis=0)
    t = np.where(denom != 0, num / denom, -np.inf)
    return 0.5 - np.arctan(t) / np.pi
    

def main(pvals_matrix, binary_matrix, fdr_trheshold=0.001):
    ones_mask = pvals_matrix == 1.
    pvals_matrix[ones_mask] = 0.5 # gets ignored in aggregation

    combined_pval = acat_equal(pvals_matrix, ones_mask)
    _, fdr, _, _ = multipletests(combined_pval, method='bonferroni')
    mcv = binary_matrix.sum(axis=0).A1
    one_pr = np.ceil(binary_matrix.shape[1] * 0.01)
    core = (mcv > one_pr) & (fdr < fdr_trheshold)
    print(core.shape, core.sum())
    return core


if __name__ == "__main__":
    samples_meta = pd.read_table(sys.argv[1]).set_index('ag_id')
    grouping_column = sys.argv[2]
    value = sys.argv[3]
    samples_meta[grouping_column] = samples_meta[grouping_column].astype(str)
    assert value in samples_meta[grouping_column].unique()
    samples = samples_meta.query(f'{grouping_column} == "{value}"').index
    anndata = read_zarr_backed(sys.argv[4])
    print('Finished reading anndata')
    mask = anndata.obs_names.isin(samples)


    pvals_matrix = np.load(sys.argv[5], mmap_mode='r')[:, mask].astype('float32')
    pvals_matrix = np.power(10, -pvals_matrix)
    # Finished reading matrix

    binary = anndata[mask, :].layers['binary'].T

    assert pvals_matrix.shape == binary.shape
    print(pvals_matrix.shape)

    fdr = float(sys.argv[6])
    core_set_mask = main(pvals_matrix, binary, fdr_trheshold=fdr)

    anndata.var[core_set_mask].reset_index()[['#chr', 'start', 'end', 'dhs_id']].to_csv(
        sys.argv[7],
        sep='\t',
        index=False,
    )

    np.save(
        sys.argv[8],
        core_set_mask
    )


    
