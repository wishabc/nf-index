import sys
import numpy as np
import json
from lowess import DataNormalize
import h5py 

def get_interpolation_for_gini(x, lowess_est, sampled):
    sort_ind = np.argsort(x[sampled])
    x_ref = np.array([x.min()] + list(x[sampled][sort_ind]) + [x.max()])
    y_ref = np.array([1] + list(lowess_est[sort_ind]) + [0])
    interpolated = np.interp(x, x_ref, y_ref)
    return interpolated


def main(normalized_matrix, binary_matrix, num_peaks, min_peaks_per_sample, meta_labels, save=None):
    unique_labels = np.unique(meta_labels)
    new_norm_matrix = np.zeros(shape=(normalized_matrix.shape[0], unique_labels.shape[0]),
                                dtype=normalized_matrix.dtype)
    for label in unique_labels:
        new_norm_matrix[:, label] = normalized_matrix[:, meta_labels == label].mean(axis=1)

    gini = np.cumsum(np.sort(new_norm_matrix, axis=1) - new_norm_matrix.min(axis=1)[:, None], axis=1)
    q = np.linspace(0, 1, new_norm_matrix.shape[1])
    gini = 2 * (q[None, :] - gini / gini[:, -1:]).mean(axis=1)
    means = np.mean(new_norm_matrix, axis=1)

    data_norm = DataNormalize(jobs=1, sample_method='raw')
    sampled_peaks_mask = data_norm.select_peaks_uniform(means, np.ones(means.shape[0], dtype=bool))
    smoothed_gini_samp = data_norm.run_lowess(gini, means, sampled_peaks_mask, frac=0.2, delta=0.001)

    smoothed_gini_final = get_interpolation_for_gini(means, smoothed_gini_samp, sampled_peaks_mask)
    gini_index = gini - smoothed_gini_final
    gini_argsort = np.argsort(gini_index)[::-1]

    top_gini_mask = gini_index > gini_index[gini_argsort[num_peaks]]
    if min_peaks_per_sample is not None and min_peaks_per_sample != 0:
        for sample_id in range(normalized_matrix.shape[1]):
            sample_binary_mask = binary_matrix[:, sample_id]
            chosen_peaks_mask = top_gini_mask * sample_binary_mask
            to_add_peaks = min_peaks_per_sample - chosen_peaks_mask.sum()
            if to_add_peaks <= 0:
                continue
            print(f"Adding peaks for {sample_id}")
            peaks_pool = sample_binary_mask * ~top_gini_mask
            if peaks_pool.sum() < to_add_peaks:
                print(f'Not enough peaks to add for {sample_id}, adding {peaks_pool.sum()}/{to_add_peaks}')
                top_gini_mask[peaks_pool] = 1
                continue
            
            to_add_argsort = peaks_pool * gini_argsort
            top_gini_mask[to_add_argsort[to_add_argsort != 0][:to_add_peaks]] = 1
    if save is not None:
        np.save(f'{save}.npy',
            normalized_matrix[
                (top_gini_mask), :
            ]
        )
    return gini_argsort, top_gini_mask



if __name__ == '__main__':
    with open(sys.argv[1], 'r') as json_file:
        params = json.load(json_file)
    with h5py.File(sys.argv[2]) as f:
        normalized_matrix = f['vsd']
        binary_matrix = f['binary']
    prefix = sys.argv[3]
    meta_labels = sys.argv[4]
    num_peaks = params['num_peaks']
    min_peaks_per_sample = params.get('min_peaks_per_sample')
    _, new_mask = main(normalized_matrix, binary_matrix, num_peaks, min_peaks_per_sample,
                    meta_labels=meta_labels, save=prefix)
    
    np.savetxt(f'{prefix}.mask.txt', new_mask, fmt="%5i")
    