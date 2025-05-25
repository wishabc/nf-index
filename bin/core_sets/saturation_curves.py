import numpy as np
from tqdm import tqdm
from genome_tools.data.anndata import read_zarr_backed
import sys


from numba import njit, prange




s_curves = {}
step_added = {}
s_curves_core = {}

n_shuffles = 100
order = np.argsort(fdr_core_set_meta['core_set_size'] * fdr_core_set_meta['n_samples'])[::-1]

for name, row in tqdm(fdr_core_set_meta.iloc[order].iterrows(), total=len(fdr_core_set_meta)):
    if s_curves_core.get(name) is not None:
        continue
    sample_mask = dnase_anndata.obs.eval(f'core_ontology_term == "{name}"').values
    mcv_mask = binary_mat[sample_mask, :].sum(axis=0) > 0
    b_mat = np.ascontiguousarray(binary_mat[sample_mask, :][:, mcv_mask])
    s_curves[name], step_added[name] = saturation_curve_dhs_add_steps(b_mat, n_shuffles)
    s_curves_core[name], _ = saturation_curve_dhs_add_steps(np.ascontiguousarray(b_mat[:, cores[row['ID']][mcv_mask]]), n_shuffles)


if __name__ == "__main__":
    anndata = read_zarr_backed(sys.argv[1])
    sparse_binary = anndata.layers['binary'].T




    

for name, row in tqdm(fdr_core_set_meta.iloc[order].iterrows(), total=len(fdr_core_set_meta)):
    if row['n_samples'] < 50:
        continue
    if meds_mcv_by_step.get(name) is not None and last_step_mcv_arrays.get(name) is not None:
        continue
    inv_sample_mask = dnase_anndata.obs.eval(f'core_ontology_term != "{name}"').values
    sample_mask = dnase_anndata.obs.eval(f'core_ontology_term == "{name}"').values
    mcv_mask = (sparse_binary[sample_mask, :].sum(axis=0) > 0).A1
    inv_b_mat = sparse_binary[inv_sample_mask, :][:, mcv_mask]
    inv_mcv = inv_b_mat.sum(axis=0).A1
    core_mask = (~cores[row['ID']]).astype(int)
    meds_mcv_by_step[name] = np.zeros(sample_mask.sum())
    q1_mcv_by_step[name] = np.zeros(sample_mask.sum())
    q3_mcv_by_step[name] = np.zeros(sample_mask.sum())
    for step in np.arange(sample_mask.sum()):
        # print(step, step_added[name].shape)
        weights = (step_added[name] == step).sum(axis=0) * core_mask[mcv_mask]
        # print('w', weights.shape, inv_mcv.shape)
        meds_mcv_by_step[name][step] = weighted_q(inv_mcv, weights)
        # print('m')
        q1_mcv_by_step[name][step] = weighted_q(inv_mcv, weights, 0.25)
        q3_mcv_by_step[name][step] = weighted_q(inv_mcv, weights, 0.75)


    last_step_mcv_arrays[name] = inv_mcv[np.where(step_added[name] == step)[1]]