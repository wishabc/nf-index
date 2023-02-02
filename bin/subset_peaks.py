import sys
import numpy as np
import json
from lowess import DataNormalize


def get_interpolation_for_gini(x, lowess_est, sampled):
    sort_ind = np.argsort(x[sampled])
    x_ref = np.array([x.min()] + list(x[sampled][sort_ind]) + [x.max()])
    y_ref = np.array([1] + list(lowess_est[sort_ind]) + [0])
    interpolated = np.interp(x, x_ref, y_ref)
    return interpolated


def main(params, normalized_matrix):
    gini = np.cumsum(np.sort(normalized_matrix, axis=1) - normalized_matrix.min(axis=1)[:, None], axis=1)
    q = np.linspace(0, 1, normalized_matrix.shape[1])
    gini = 2 * (q[None, :] - gini / gini[:, -1:]).mean(axis=1)
    means = np.mean(normalized_matrix, axis=1)

    data_norm = DataNormalize(jobs=1)
    sampled_peaks_mask = data_norm.select_peaks_uniform(means, singletons_mask,
                                                       sample_method='raw')
    smoothed_gini_samp = data_norm.run_lowess(gini, means, sampled_peaks_mask,
                                                frac=0.2, delta=0.001)

    smoothed_gini_final = get_interpolation_for_gini(means, smoothed_gini_samp, sampled_peaks_mask)
    gini_index = gini - smoothed_gini_final

    num_peaks = params['num_peaks']
    np.save(out_path,
        normalized_matrix[
            (gini_index > gini_index[np.argsort(gini_index)[::-1][num_peaks]]), :
        ]
    )


if __name__ == '__main__':
    with open(sys.argv[1], 'r') as json_file:
        params = json.load(json_file)
    
    normalized_matrix = np.load(sys.argv[2])
    singletons_mask = np.loadtxt(sys.argv[3]).astype(bool)
    out_path = sys.argv[4]
    assert normalized_matrix.shape[0] == singletons_mask.shape[0]
    normalized_matrix = normalized_matrix[singletons_mask, :]
    main(params, normalized_matrix)
    