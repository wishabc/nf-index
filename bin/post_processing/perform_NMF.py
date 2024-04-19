import numpy as np
import pandas as pd
import argparse

from sklearn.decomposition import NMF
from weighted_NMF import NMF as weighted_NMF


def perform_NMF(X, weights=None, n_components=16, model=None):
    if model is None:
        params = dict(n_components=n_components,
                        solver='mu', beta_loss='frobenius',
                        random_state=0, init="nndsvda",
                        max_iter=1000, tol=1e-4, verbose=True)
        if weights is not None:
            model = weighted_NMF(**params)
            W = model.fit_transform(X.T, weights=weights)
        else:
            model = NMF(**params)
            W = model.fit_transform(X.T)
    else:
        if weights is not None:
            W = model.fit_transform(X.T, weights=weights)
        else:
            W = model.transform(X.T)
    H = model.components_ # components x peaks
    return W, H, model

def project_samples(data, model):
    # data: peaks x samples 
    # NMF: samples x peaks = samples x components * components x peaks
    return model.transform(data.T) # samples x components 

def project_peaks(data, model, W):
    # data: peaks x samples
    # NMF: peaks x samples = peaks x components * components x samples
    return model._fit_transform(data, H=W, update_H=False) # components x peaks


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Matrix normalization using lowess')
    parser.add_argument('matrix', help='Path to matrix to run NMF on')
    parser.add_argument('prefix', help='Prefix for the output file')
    parser.add_argument('n_components', help='Number of components to use in NMF', type=int)
    parser.add_argument('--samples_mask', help='Mask of used samples, numpy array', default=None)
    parser.add_argument('--peaks_mask', help='Mask of used samples, numpy array', default=None)
    parser.add_argument('--samples_weights', help='Path to samples weights (for weighted NMF)', default=None)
    args = parser.parse_args()

    mat = np.load(args.matrix)
    
    if args.samples_mask is not None:
        samples_m = np.load(args.samples_mask)
    else:
        samples_m = np.ones(mat.shape[1], dtype=bool)
    
    if args.peaks_mask is not None:
        peaks_m = np.loadtxt(args.peaks_mask).astype(bool)
    else:
        peaks_m = np.ones(mat.shape[0], dtype=bool)
    
    weights_vector = None
    if args.samples_weights:
        weights_df = pd.read_table(args.samples_weights)
        weights_vector = weights_df.set_index("id").to_numpy().squeeze()
        if args.samples_mask is not None:
            weights_vector = weights_vector[samples_m]

    mat = mat[peaks_m, :].astype(float)
    matrix = mat[:, samples_m]
    non_zero_rows = matrix.sum(axis=1) > 0

    matrix = matrix[non_zero_rows, :]

    print('Fitting NMF model')
    W_np, H_np, model = perform_NMF(X=matrix, weights=weights_vector, n_components=args.n_components)
    if args.samples_mask is not None:
        print('Projecting samples')
        W_np = project_samples(mat[non_zero_rows, :], model)
        if non_zero_rows.sum() != mat.shape[0]:
            H_np = project_peaks(mat, model, W_np)
    
    print('Saving results')
    np.save(f'{args.prefix}.W.npy', W_np) # samples x components
    np.save(f'{args.prefix}.H.npy', H_np) # components x peaks
    np.savetxt(f'{args.prefix}.non_zero_peaks_mask.txt', non_zero_rows, fmt="%d")
    np.savetxt(f'{args.prefix}.samples_mask.txt', samples_m, fmt="%d")