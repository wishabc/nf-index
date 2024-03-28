import sys
import numpy as np
import json
from ..lowess import DataNormalize
import h5py 

def get_interpolation_for_gini(x, lowess_est, sampled):
    sort_ind = np.argsort(x[sampled])
    x_ref = np.array([x.min()] + list(x[sampled][sort_ind]) + [x.max()])
    y_ref = np.array([1] + list(lowess_est[sort_ind]) + [0])
    interpolated = np.interp(x, x_ref, y_ref)
    return interpolated


def average_matrices(normalized_matrix, binary_matrix, meta_labels, reprod=1):
    unique_labels = np.unique(meta_labels)
    new_norm_matrix = np.zeros(shape=(normalized_matrix.shape[0], unique_labels.shape[0]),
                                dtype=normalized_matrix.dtype)
    new_binary_matrix = np.zeros(shape=(binary_matrix.shape[0], unique_labels.shape[0]),
                                dtype=bool)
    for label in unique_labels:
        new_norm_matrix[:, label] = normalized_matrix[:, meta_labels == label].mean(axis=1)
        new_binary_matrix[:, label] = binary_matrix[:, meta_labels == label].sum(axis=1) >= np.minimum(reprod, (meta_labels == label).sum())

    return new_norm_matrix, new_binary_matrix


def get_gini_index_for_peaks(new_norm_matrix, num_peaks):
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

    return means, gini, smoothed_gini_final, gini_argsort, top_gini_mask


def add_peaks(binary_matrix, gini_argsort, num_peaks=0, min_peaks_per_sample=0):
    sorted_binary = binary_matrix[gini_argsort, :]
    sorted_binary = np.where(sorted_binary.cumsum(axis=0) > min_peaks_per_sample, False, sorted_binary)
    sorted_binary[:num_peaks, :] = True

    inv_argsort = np.zeros_like(gini_argsort)
    inv_argsort[gini_argsort] = np.arange(len(gini_argsort))

    return sorted_binary.any(axis=1)[inv_argsort]


def main(normalized_matrix, binary_matrix, num_peaks, min_peaks_per_sample, meta_labels,
    meta_labels_in_gini=True, meta_labels_in_add=False, save=None):
    new_norm_matrix, new_binary_matrix = average_matrices(normalized_matrix, binary_matrix, meta_labels)

    if meta_labels_in_gini:
        means, gini, smoothed_gini_final, gini_argsort, top_gini_mask = get_gini_index_for_peaks(new_norm_matrix, num_peaks)
    else:
        means, gini, smoothed_gini_final, gini_argsort, top_gini_mask = get_gini_index_for_peaks(normalized_matrix, num_peaks)

    if meta_labels_in_add:
        top_gini_mask = add_peaks(new_binary_matrix, gini_argsort, num_peaks, min_peaks_per_sample)
    else:
        top_gini_mask = add_peaks(binary_matrix, gini_argsort, num_peaks, min_peaks_per_sample)

    if save is not None:
        np.save(f'{save}.npy',
            normalized_matrix[
                (top_gini_mask), :
            ]
        )
    return gini_argsort, top_gini_mask, means, gini, smoothed_gini_final



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
    _, new_mask, *_ = main(normalized_matrix, binary_matrix, num_peaks, min_peaks_per_sample,
                    meta_labels=meta_labels, save=prefix)
    
    np.savetxt(f'{prefix}.mask.txt', new_mask, fmt="%5i")
    