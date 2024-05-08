import numpy as np
import pandas as pd
import argparse

from sklearn.decomposition import NMF
from weighted_NMF import NMF as weighted_NMF


def perform_NMF(X, W_weights=None, H_weights=None, n_components=16):
    params = dict(n_components=n_components,
                    solver='mu', beta_loss='frobenius',
                    random_state=0, init="nndsvda",
                    max_iter=1000, tol=1e-4, verbose=True)
    if W_weights is not None:
        assert H_weights is not None
        model = weighted_NMF(**params)
        W = model.fit_transform(X.T, W_weights=W_weights[:, None], H_weights=H_weights[None, :])
    else:
        model = NMF(**params)
        W = model.fit_transform(X.T)

    H = model.components_ # components x peaks
    return W, H, model

def project_samples(data, model, H):
    # data: peaks x samples
    # H: components x peaks
    # NMF: samples x peaks = samples x components * components x peaks
    W, *_ = model._fit_transform(data.T, H=H, update_H=False)
    return W # samples x components 

def project_peaks(data, model, W, W_weights=None, H_weights=None):
    # data: peaks x samples
    # W: samples x components
    # NMF: peaks x samples = peaks x components * components x samples
    params = dict(X=data, H=W.T, update_H=False)
    if W_weights is not None:
        assert H_weights is not None
        projected_peaks, _, _ = model._fit_transform(
            **params,
            H_weights=W_weights[None, :],
            W_weights=H_weights[:, None]
        ) # components x peaks
    else:
        projected_peaks, _, _ = model._fit_transform(**params)
    return projected_peaks.T

def get_nonzero_mask(matrix):
    return matrix.sum(axis=1) > 0

def read_weights(weights_path, shape, ext=None):
    weights_vector = np.ones(shape, dtype=float)
    if weights_path:
        if ext == 'npy':
            weights_vector = np.load(weights_path)  
        else:
            weights_df = pd.read_table(weights_path)
            weights_vector = weights_df.set_index("ag_id").to_numpy().squeeze()
    
    return weights_vector / weights_vector.sum() * weights_vector.shape[0]

def read_args(args):
    print('Reading matrix')
    mat = np.load(args.matrix).astype(float)
    
    if args.samples_mask is not None:
        samples_m = np.load(args.samples_mask)
    else:
        samples_m = np.ones(mat.shape[1], dtype=bool)
    
    if args.peaks_mask is not None:
        peaks_m = np.load(args.peaks_mask).astype(bool)
    else:
        peaks_m = np.ones(mat.shape[0], dtype=bool)
    
    
    if args.samples_weights or args.peaks_weights:
        W_weights_vector = read_weights(args.samples_weights, mat.shape[1])
        H_weights_vector = read_weights(args.peaks_weights, mat.shape[0])
    else:
        H_weights_vector = W_weights_vector = None

    return mat, samples_m, peaks_m, W_weights_vector, H_weights_vector


def main(mat, samples_m, peaks_m, W_weights, H_weights):
    samples_masked_matrix = mat[:, samples_m]
    non_zero_rows = samples_masked_matrix.sum(axis=1) > 0

    peaks_mask = peaks_m & non_zero_rows

    matrix_samples_peaks_slice = mat[peaks_mask, :][:, samples_m]

    print('Fitting NMF model')
    if W_weights is not None or H_weights is not None:
        print('Using weighted NMF')
        W_weights_slice = W_weights[samples_m]
        H_weights_slice = H_weights[peaks_m]

    else:
        W_weights_slice = H_weights_slice = None

    W_np, H_np, model = perform_NMF(
        X=matrix_samples_peaks_slice,
        W_weights=W_weights_slice,
        H_weights=H_weights_slice,
        n_components=args.n_components
    )
    if samples_m.shape[0] > samples_m.sum():
        assert W_weights is None
        print('Projecting samples')
        if args.project_masked_peaks and (non_zero_rows.sum() != mat[:, samples_m].shape[0]):
            H_np = project_peaks(mat[:, samples_m][non_zero_rows, :], model, W_np)
            W_np = project_samples(mat[non_zero_rows, :], model, H_np)
        else:
            W_np = project_samples(mat[peaks_mask, :], model, H_np)

    if peaks_mask.shape[0] > peaks_mask.sum():
        H_np = project_peaks(mat, model, W_np, W_weights, H_weights)

    return W_np, H_np, peaks_mask

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Matrix normalization using lowess')
    parser.add_argument('matrix', help='Path to matrix to run NMF on')
    parser.add_argument('prefix', help='Prefix for the output file')
    parser.add_argument('n_components', help='Number of components to use in NMF', type=int)
    parser.add_argument('--samples_mask', help='Mask of used samples, numpy array', default=None)
    parser.add_argument('--project_masked_peaks', action='store_true', help='Project peaks for all masked samples', default=False)
    parser.add_argument('--peaks_mask', help='Mask of used samples, numpy array', default=None)
    parser.add_argument('--samples_weights', help='Path to samples weights (for weighted NMF)', default=None)
    parser.add_argument('--peaks_weights', help='Path to peaks weights (for weighted NMF)', default=None)
    args = parser.parse_args()

    mat, samples_m, peaks_m, W_weights_vector, H_weights_vector = read_args(args)
    W_np, H_np, peaks_mask = main(mat, samples_m, peaks_m, W_weights_vector, H_weights_vector)
    
    print('Saving results')
    np.save(f'{args.prefix}.W.npy', W_np) # samples x components
    np.save(f'{args.prefix}.H.npy', H_np.T) # components x peaks
    np.savetxt(f'{args.prefix}.non_zero_peaks_mask.txt', peaks_mask, fmt="%d")
    np.savetxt(f'{args.prefix}.samples_mask.txt', samples_m, fmt="%d")
    #np.savetxt(f'{args.prefix}.weights_vector.txt', weights_vector)