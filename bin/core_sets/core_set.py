import argparse
import numpy as np
import pandas as pd
from genome_tools.data.anndata import read_zarr_backed
from statsmodels.stats.multitest import multipletests

from scipy.stats import cauchy, combine_pvalues, norm
from numba import njit, prange
from tqdm import tqdm


def weighted_q(values, weights, q=0.5):
    """
    Compute the weighted median of values with corresponding weights.
    For integer weights this is equivalent to repeating each value weight times
    and taking the ordinary median.
    """
    a = np.asarray(values)
    w = np.asarray(weights)
    # sort both arrays by value
    idx = np.argsort(a)
    a_sorted, w_sorted = a[idx], w[idx]
    # cumulative weight and cutoff
    cumw = np.cumsum(w_sorted)
    cutoff = w_sorted.sum() * q
    # first value where cum weight ≥ q total weight
    return a_sorted[cumw >= cutoff][0]

def stouffer_equal_weights(p, ones_mask, denom):
    """
    Equal-weight Stouffer’s method on p-value matrix `p`, masked by `ones_mask`.
    Returns combined p-values per column.
    """
    z = norm.isf(p)
    z = np.where(ones_mask, 0.0, z)

    z_comb = np.sum(z, axis=0) / np.sqrt(denom)
    return norm.sf(z_comb).astype(p.dtype)


def acat_equal_weights(p, ones_mask, denom):
    q = cauchy.ppf(p)
    q = np.where(ones_mask, 0.0, q)
    t = np.sum(q, axis=0) / denom
    return cauchy.cdf(t).astype(p.dtype)


@njit(parallel=True)
def saturation_curve_dhs_add_steps(binary_mat, num_shuffles):
    N, M = binary_mat.shape
    cum_sum = np.empty((num_shuffles, N), dtype=np.int64)
    dhs_add_step = -np.ones((num_shuffles, M), dtype=np.int64)  # -1 = never contributed

    for s in prange(num_shuffles):
        perm = np.arange(N, dtype=np.int64)
        for j in range(N - 1, 0, -1):
            k = np.random.randint(0, j + 1)
            perm[j], perm[k] = perm[k], perm[j]

        current = np.zeros(M, dtype=np.bool_)
        count = 0

        for step in range(N):
            row = binary_mat[perm[step]]
            for j in range(M):
                if row[j] and not current[j]:
                    current[j] = True
                    count += 1
                    dhs_add_step[s, j] = step
            cum_sum[s, step] = count

    return cum_sum, dhs_add_step


def calc_mcv(binary_sparse):
    return binary_sparse.sum(axis=0).A1


def saturation_curve(category_binary, mcv_mask, core_mask, n_shuffles=100):
    b_mat = np.ascontiguousarray(category_binary[:, mcv_mask].todense().A)
    core_b_mat = np.ascontiguousarray(category_binary[:, core_mask].todense().A)
    print('Saturation curve', b_mat.shape, core_b_mat.shape, n_shuffles)
    s_curve, step_added = saturation_curve_dhs_add_steps(b_mat, n_shuffles)
    s_curve_core, _ = saturation_curve_dhs_add_steps(core_b_mat, n_shuffles)
    print('Finished saturation curve', b_mat.shape, core_b_mat.shape, n_shuffles)
    return s_curve, s_curve_core, step_added


def per_step_stats(category_binary, inv_mcv, step_added, mcv_mask, core_mask):
    n_samples = category_binary.shape[0]
    mcv_by_step_stats = np.zeros((n_samples, 3))  # median, q1, q3

    for step in tqdm(np.arange(n_samples)):
        # print(step, step_added[name].shape)
        weights = (step_added == step).sum(axis=0) * core_mask[mcv_mask]
        # print('w', weights.shape, inv_mcv.shape)
        median = weighted_q(inv_mcv, weights)
        q1 = weighted_q(inv_mcv, weights, 0.25)
        q3 = weighted_q(inv_mcv, weights, 0.75)
        mcv_by_step_stats[step, :] = median, q1, q3

    return mcv_by_step_stats


def combine_pvals(pvals_matrix, combine_pval_method):
    ones_mask = pvals_matrix == 1.
    pvals_matrix = np.where(ones_mask, 0.5, pvals_matrix)  # gets ignored in aggregation
    denom = (~ones_mask).sum(axis=0)
    if combine_pval_method == 'cauchy':
        combined_pval = acat_equal_weights(pvals_matrix, ones_mask, denom)
    elif combine_pval_method == 'stouffer':
         combined_pval = stouffer_equal_weights(pvals_matrix, ones_mask, denom)
    else:
        raise ValueError(f'Unknown method for p-value combination: {combine_pval_method}')
    combined_pval[np.isnan(combined_pval)] = 1.0
    return combined_pval

def get_core_set(pvals_matrix, binary_matrix, combine_pval_method, fdr_threshold=0.001):
    combined_pval = combine_pvals(pvals_matrix, combine_pval_method=combine_pval_method)
    fdr = multipletests(combined_pval, method='bonferroni')[1]
    print(fdr.min())
    mcv = binary_matrix.astype(int).sum(axis=0).A1
    one_pr = np.ceil(binary_matrix.shape[0] * 0.01)
    print(one_pr, (mcv > one_pr).sum())

    core = (mcv >= one_pr) & (fdr <= fdr_threshold)
    print(core.shape, core.sum())
    return core


def main(pvals_matrix, binary, inverse_binary, fdr_threshold, combine_pval_method='cauchy'):
    category_binary = binary
    mcv_mask = calc_mcv(category_binary) > 0
    inv_mcv = calc_mcv(inverse_binary[:, mcv_mask])
    core_mask = get_core_set(
        pvals_matrix,
        category_binary,
        fdr_threshold=fdr_threshold,
        combine_pval_method=combine_pval_method
    )
    s_curve, s_curve_core, step_added = saturation_curve(
        category_binary,
        mcv_mask,
        core_mask
    )
    mcv_by_step_stats = per_step_stats(
        category_binary,
        inv_mcv,
        step_added,
        mcv_mask,
        core_mask
    )

    return core_mask, s_curve, s_curve_core, step_added, mcv_by_step_stats


def parse_args():
    parser = argparse.ArgumentParser(description="Compute core sets from p-values and binary matrix along with saturation curves.")

    parser.add_argument("grouping_column", type=str, help="Column in samples_meta to use for grouping.")
    parser.add_argument("output_prefix", type=str, help="Prefix for all output files.")
    parser.add_argument("value", type=str, help="Value in grouping_column to select samples.")
    parser.add_argument("anndata_path", type=str, help="Path to Zarr-backed AnnData object.")
    parser.add_argument("--fdr", type=float, help="FDR threshold for core set selection.")
    parser.add_argument("--method", type=str, choices=("cauchy", "stouffer"), help="Method for p-value combination.")

    return parser.parse_args()


def load_and_validate_inputs(grouping_column, value, anndata_path, pvals_path):
    anndata = read_zarr_backed(anndata_path)
    samples_meta = anndata.obs
    samples_meta[grouping_column] = samples_meta[grouping_column].astype(str)

    if value not in samples_meta[grouping_column].unique():
        raise ValueError(f"Value '{value}' not found in column '{grouping_column}'.")

    samples_mask = samples_meta.eval(f'{grouping_column} == "{value}"').values
    inverse_binary = anndata[~samples_mask, :].layers['binary']
    anndata = anndata[samples_mask, :]
    print("Finished reading anndata.")

    return anndata, inverse_binary


if __name__ == "__main__":
    args = parse_args()

    anndata, inverse_binary = load_and_validate_inputs(
        args.grouping_column,
        args.value,
        args.anndata_path,
    )
    binary = anndata.layers["binary"]
    pvals_matrix = anndata.layers['neglog10_pvals'].astype(np.float64)
    pvals_matrix = np.power(10, -pvals_matrix)  # Convert -log10(p) to p

    core_set_mask, s_curve, s_curve_core, step_added, mcv_by_step_stats = main(
        pvals_matrix,
        binary,
        inverse_binary,
        fdr_threshold=args.fdr,
        combine_pval_method=args.method
    )

    # Output files
    prefix = args.output_prefix

    anndata.var[core_set_mask].reset_index()[['#chr', 'start', 'end', 'dhs_id']].to_csv(
        f"{prefix}.core_set.bed", sep='\t', index=False
    )
    np.save(f"{prefix}.core_set.npy", core_set_mask)
    np.save(f"{prefix}.saturation_curve.npy", s_curve)
    np.save(f"{prefix}.saturation_curve_core.npy", s_curve_core)
    np.save(f"{prefix}.step_added.npy", step_added)
    np.save(f"{prefix}.mcv_by_step_stats.npy", mcv_by_step_stats)
