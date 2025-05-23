import sys
import numpy as np
import pandas as pd
from genome_tools.data.anndata import read_zarr_backed
from statsmodels.stats.multitest import multipletests
from scipy.stats import cauchy

def acat_equal(p, ones_mask):
    # z = norm.isf(p)
    # z_sum = np.sum(z, axis=0)
    # z_comb = z_sum / np.sqrt(z.shape[0])
    # return norm.sf(z_comb)
    denom = (~ones_mask).sum(axis=0)
    q = cauchy.ppf(p)
    t = np.sum(q, axis=0) / denom
    return cauchy.cdf(t).astype(p.dtype)


def main(pvals_matrix, binary_matrix, fdr_threshold=0.001):
    ones_mask = pvals_matrix == 1.
    pvals_matrix = np.where(ones_mask, 0.5, pvals_matrix)  # gets ignored in aggregation
    combined_pval = acat_equal(pvals_matrix, ones_mask)
    combined_pval[np.isnan(combined_pval)] = 1.0
    fdr = multipletests(combined_pval, method='bonferroni')[1]
    print(fdr.min())
    mcv = binary_matrix.astype(int).sum(axis=0).A1
    one_pr = np.ceil(binary_matrix.shape[0] * 0.01)
    print(one_pr, (mcv > one_pr).sum())

    core = (mcv >= one_pr) & (fdr <= fdr_threshold)
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


    pvals_matrix = np.load(sys.argv[5], mmap_mode='r')[:, mask].astype('float32').T
    pvals_matrix = np.power(10, -pvals_matrix)
    # Finished reading matrix

    binary = anndata[mask, :].layers['binary']

    assert pvals_matrix.shape == binary.shape
    print(pvals_matrix.shape)

    fdr_tr = float(sys.argv[6])
    core_set_mask = main(pvals_matrix, binary, fdr_threshold=fdr_tr)

    anndata.var[core_set_mask].reset_index()[['#chr', 'start', 'end', 'dhs_id']].to_csv(
        sys.argv[7],
        sep='\t',
        index=False,
    )

    np.save(
        sys.argv[8],
        core_set_mask
    )


    
