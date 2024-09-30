import numpy as np
import pandas as pd
import argparse
import json
from sklearn.decomposition import NMF
from weighted_NMF import WeightedNMF
from collections import namedtuple
from nmf_tools.data.anndata_loader import read_zarr_backed
import os

NMFInputData = namedtuple('NMFInputData', ['matrix', 'samples_mask', 'peaks_mask', 'samples_weights', 'peaks_weights', 'samples_metadata', 'dhs_metadata'])


def initialize_model(n_components, extra_params, is_weighted=False):
    params = dict(n_components=n_components,
                    solver='mu', beta_loss='frobenius',
                    random_state=0, init="nndsvda",
                    max_iter=1000, tol=1e-4, alpha_W=0.0, l1_ratio=1.0,
                    verbose=True)
    if extra_params is not None and extra_params:
        overwritten = {
            key: f'New: {extra_params[key]}. Old: {params[key]}' for key in extra_params
              if key in params and extra_params[key] != params[key]
        }
        params.update(extra_params)
        print("Overwritten params:", overwritten)
    
    return WeightedNMF(**params) if is_weighted else NMF(**params)


def run_NMF(model: NMF, X, W_weights=None, H_weights=None):
    assert W_weights is None == H_weights is None
    if W_weights is not None:
        assert len(W_weights.shape) == 1 and len(H_weights.shape) == 1, 'Weights should be 1D arrays'
        W = model.fit_transform(X.T, W_weights=W_weights[:, None], H_weights=H_weights[None, :])
    else:
        W = model.fit_transform(X.T)

    H = model.components_ # components x peaks
    return W, H, model


def project_samples(data, model, H, W_weights=None, H_weights=None):
    # data: peaks x samples
    # H: components x peaks
    # NMF: samples x peaks = samples x components * components x peaks
    params = dict(X=data.T, H=H, update_H=False)
    if W_weights is not None:
        assert H_weights is not None
        params = {**params,
            'H_weights': H_weights[None, :],
            'W_weights': W_weights[:, None]
        }
    W, *_ = model._fit_transform(**params)
    return W # samples x components 


def project_peaks(data, model, W, W_weights=None, H_weights=None):
    # data: peaks x samples
    # W: samples x components
    # NMF: peaks x samples = peaks x components * components x samples
    params = dict(X=data, H=W.T, update_H=False)
    if W_weights is not None:
        assert H_weights is not None
        params = {**params,
            'H_weights': H_weights[None, :],
            'W_weights': W_weights[:, None]
        }
    projected_peaks, _, _ = model._fit_transform(**params)
    return projected_peaks.T


##################
### Parse args ###
##################

def parse_nmf_args(args) -> NMFInputData:
    if args.from_anndata is not None:
        return parse_args_anndata(args)
    return parse_args_matrix(args)


def parse_args_anndata(args):
    anndata = read_zarr_backed(args.from_anndata)
    anndata = anndata[anndata.obs['final_qc_passing_sample'], anndata.var['final_qc_passing_dhs']]
    matrix = anndata.layers['binary'].T.todense()
    return parse_optional_args(args, matrix, anndata.obs, anndata.var)


def parse_args_matrix(args):
    print('Reading matrix')
    mat = np.load(args.matrix).astype(np.float32)

    sample_names = np.loadtxt(args.sample_names, dtype=str)

    try:
        dhs_meta = pd.read_table(args.dhs_meta, header=None, usecols=np.arange(4), names=['chr', 'start', 'end', 'dhs_id']).set_index('dhs_id')
        assert dhs_meta.shape[0] == mat.shape[0]
    except AssertionError:
        dhs_meta = pd.read_table(args.dhs_meta)
        assert dhs_meta.shape[0] == mat.shape[0]

    if hasattr(args, 'samples_metadata'):
        samples_metadata = pd.read_table(args.samples_metadata)
        # Reformat metadata to handle DNase columns
        id_col = 'id' if 'id' in samples_metadata.columns else 'ag_id'
        assert id_col in samples_metadata.columns, f'No id or ag_id column found in metadata. Available columns: {samples_metadata.columns}'
    
        samples_metadata = samples_metadata.set_index(id_col).loc[sample_names]
    else:
        samples_metadata = pd.DataFrame(index=sample_names)
    
    if hasattr(args, 'dhs_annotations') \
        and (args.dhs_annotation is not None) \
        and os.path.exists(args.dhs_annotations) \
        and ('dist_tss' not in dhs_meta.columns):
        dhs_annotations = pd.read_table(args.dhs_annotations).set_index('dhs_id')
        dhs_meta['dist_tss'] = dhs_annotations['dist_tss']
    
    return parse_optional_args(args, mat, dhs_metadata=dhs_meta, samples_metadata=samples_metadata)


def parse_optional_args(args, matrix, samples_metadata, dhs_metadata):
    if args.samples_mask is not None:
        samples_m = np.load(args.samples_mask)
    else:
        samples_m = np.ones(matrix.shape[1], dtype=bool)
    
    if args.peaks_mask is not None:
        peaks_m = np.load(args.peaks_mask).astype(bool)
    else:
        peaks_m = np.ones(matrix.shape[0], dtype=bool)
    
    if args.samples_weights or args.peaks_weights:
        W_weights_vector = read_weights(args.samples_weights, matrix.shape[1], samples_metadata.index)
        H_weights_vector = read_weights(args.peaks_weights, matrix.shape[0], dhs_metadata.index)
    else:
        H_weights_vector = W_weights_vector = None
    
    return NMFInputData(
        matrix=matrix,
        samples_mask=samples_m,
        peaks_mask=peaks_m,
        samples_weights=W_weights_vector,
        peaks_weights=H_weights_vector,
        dhs_metadata=dhs_metadata,
        samples_metadata=samples_metadata
    )


def read_weights(weights_path, shape, sample_names=None, ext=None):
    weights_vector = np.ones(shape, dtype=float)
    if weights_path:
        if ext == 'npy':
            weights_vector = np.load(weights_path)  
        else:
            weights_df = pd.read_table(weights_path).set_index('id')
            if sample_names is not None:
                weights_df = weights_df.loc[sample_names]
            weights_vector = weights_df['weight'].to_numpy()
    
    return weights_vector / weights_vector.sum() * weights_vector.shape[0]


# Internal function to run NMF decomposition as part of pipeline
def main(nmf_input_data: NMFInputData, **extra_params):
    mat, samples_m, peaks_m, W_weights, H_weights, *_ = nmf_input_data

    samples_masked_matrix = mat[:, samples_m]
    non_zero_rows = samples_masked_matrix.sum(axis=1) > 0

    peaks_mask = peaks_m & non_zero_rows

    matrix_samples_peaks_slice = mat[peaks_mask, :][:, samples_m]

    print('Fitting NMF model')
    if W_weights is not None or H_weights is not None:
        print('Using weighted NMF')
        W_weights_slice = W_weights[samples_m]
        H_weights_slice = H_weights[peaks_m]
        model = initialize_model(n_components=args.n_components, extra_params=extra_params, is_weighted=True)
    else:
        W_weights_slice = H_weights_slice = None
        model = initialize_model(n_components=args.n_components, extra_params=extra_params, is_weighted=False)

    W_np, H_np, model = run_NMF(
        model=model,
        X=matrix_samples_peaks_slice,
        W_weights=W_weights_slice,
        H_weights=H_weights_slice,
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
    parser = argparse.ArgumentParser('Run NMF decomposition on a DHS x samples matrix')
    parser.add_argument('n_components', help='Number of components to use in NMF', type=int)
    parser.add_argument('prefix', help='Prefix for the output file')

    # Reading from npy matrix arguments
    parser.add_argument('--matrix', help='Path to matrix to run NMF on')
    parser.add_argument('--sample_names', help='Path to file with sample names')
    parser.add_argument('--dhs_meta', help='Path to DHS metadata file')

    # Reading from anndata argument
    parser.add_argument('--from_anndata', help='Path to AnnData file. If provided, ignore matrix, sample_names and dhs_meta fields', default=None)

    # Optional arguments for both modes
    parser.add_argument('--samples_mask', help='Mask of used samples, numpy array', default=None)
    parser.add_argument('--peaks_mask', help='Mask of used peaks, numpy array', default=None)
    parser.add_argument('--project_masked_peaks', action='store_true', help='Project peaks for all masked samples', default=False)
    parser.add_argument('--samples_weights', help='Path to samples weights (for weighted NMF)', default=None)
    parser.add_argument('--peaks_weights', help='Path to peaks weights (for weighted NMF)', default=None)
    parser.add_argument('--extra_params', help='Path to json with NMF params (overwrites default settings)', default=None)

    args = parser.parse_args()

    nmf_data = parse_nmf_args(args)
    if args.extra_params is not None:
        with open(args.extra_params) as f:
            extra_params = json.load(f)
    else:
        extra_params = {}

    W_np, H_np, peaks_mask = main(nmf_data, **extra_params)
    
    print('Saving results')
    np.save(f'{args.prefix}.W.npy', W_np) # samples x components
    np.save(f'{args.prefix}.H.npy', H_np.T) # components x peaks
    np.savetxt(f'{args.prefix}.non_zero_peaks_mask.txt', peaks_mask, fmt="%d")
    np.savetxt(f'{args.prefix}.samples_mask.txt', nmf_data.samples_mask, fmt="%d")
    #np.savetxt(f'{args.prefix}.weights_vector.txt', weights_vector)