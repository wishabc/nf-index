import numpy as np
import pandas as pd
import argparse
import json
from sklearn.decomposition import NMF
from weighted_NMF import WeightedNMF
from genome_tools.data.anndata import read_zarr_backed
import os
import scipy.sparse as sp
import dataclasses


dtype = np.float64


@dataclasses.dataclass
class NMFInputData:
    matrix: np.ndarray
    samples_mask: np.ndarray
    peaks_mask: np.ndarray
    samples_weights: np.ndarray
    peaks_weights: np.ndarray
    samples_metadata: pd.DataFrame
    dhs_metadata: pd.DataFrame


def initialize_model(n_components, extra_params: dict=None, is_weighted=False):
    params = dict(n_components=n_components,
                    solver='mu', beta_loss='frobenius',
                    random_state=0, init="nndsvda",
                    max_iter=1000, tol=1e-4, alpha_W=0.0, l1_ratio=1.0,
                    verbose=True)
    if extra_params:
        overwritten = {
            key: f'New: {extra_params[key]}. Old: {params[key]}' for key in extra_params
              if key in params and extra_params[key] != params[key]
        }
        params.update(extra_params)
        if overwritten:
            print("Overwritten params:", overwritten, flush=True)
    
    return WeightedNMF(**params) if is_weighted else NMF(**params)

def data_to_sparse(X: np.ndarray) -> sp.csr_matrix:
    if sp.issparse(X):
        return X.astype(dtype)
    return sp.coo_matrix(X.T).tocsr().astype(dtype)

def run_NMF(model: NMF, X, W_weights: np.ndarray=None, H_weights: np.ndarray=None):
    assert (W_weights is None) == (H_weights is None), 'Both or neither weights should be provided'
    X = data_to_sparse(X)
    if W_weights is not None:
        assert W_weights.ndim == 1 and H_weights.ndim == 1, 'Weights are expected to be 1D arrays'
        print(X.shape, W_weights.shape, H_weights.shape)
        W = model.fit_transform(X, W_weights=W_weights[:, None], H_weights=H_weights[None, :])
    else:
        W = model.fit_transform(X)

    H = model.components_ # components x peaks
    W = W.astype(X.dtype)
    H = H.astype(X.dtype)
    return W, H, model


def project_samples(model, X, H, W_weights=None, H_weights=None):
    # X: peaks x samples
    # H: components x peaks
    # NMF: samples x peaks = samples x components * components x peaks
    X = data_to_sparse(X) # transposes
    params = dict(X=X, H=H, update_H=False)
    if W_weights is not None:
        assert H_weights is not None
        params.update({
            'H_weights': H_weights[None, :],
            'W_weights': W_weights[:, None]
        })
    W, *_ = model._fit_transform(**params)
    return W # samples x components 


def project_peaks(model, X, W, W_weights=None, H_weights=None):
    # X: peaks x samples
    # W: samples x components
    # NMF: peaks x samples = peaks x components * components x samples
    X = data_to_sparse(X.T) # transposes
    params = dict(X=X, H=W.T, update_H=False)
    if W_weights is not None:
        assert H_weights is not None
        params = {**params,
            'W_weights': H_weights[:, None],
            'H_weights': W_weights[None, :]
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
    print('Reading AnnData')
    adata = read_zarr_backed(args.from_anndata)[:, :]  # to avoid csr dataset

    # ### Debug part start ###
    # indices = np.arange(adata.shape[1])
    # np.random.seed(42)
    # selected_indices = indices[np.random.rand(adata.shape[1]) < 0.05]
    # print(f'adata shape: {adata.shape[1]}')
    # selected_indices = np.sort(selected_indices)
    # np.save('/home/sabramov/tmp/selected_indices.npy', selected_indices)
    # adata = adata[:, selected_indices]
    # ### Debug part end ###

    matrix = adata.layers['binary'].T.toarray()

    if args.samples_mask is None:
        args.samples_mask = mask_from_metadata(adata.obs, args.samples_mask_column)
    
    if args.peaks_mask is None:
        args.peaks_mask = mask_from_metadata(adata.var, args.dhs_mask_column)

    result = parse_optional_args(
        args,
        matrix,
        samples_metadata=adata.obs.copy(),
        dhs_metadata=adata.var.copy()
    )
    print('Finished reading AnnData')
    return result


def parse_args_matrix(args):
    print('Reading matrix')
    mat = np.load(args.matrix)

    sample_names = np.loadtxt(args.sample_names, dtype=str)

    try:
        dhs_meta = pd.read_table(
            args.dhs_meta,
            header=None,
            usecols=np.arange(4),
            names=['chr', 'start', 'end', 'dhs_id']
        ).set_index('dhs_id')
        assert dhs_meta.shape[0] == mat.shape[0]
    except AssertionError:
        dhs_meta = pd.read_table(args.dhs_meta)
        assert dhs_meta.shape[0] == mat.shape[0]

    if args.samples_metadata is not None:
        samples_metadata = pd.read_table(args.samples_metadata)
        # Reformat metadata to handle DNase columns
        id_col = 'id' if 'id' in samples_metadata.columns else 'ag_id'
        assert id_col in samples_metadata.columns, f'No id or ag_id column found in metadata. Available columns: {samples_metadata.columns}'
    
        samples_metadata = samples_metadata.set_index(id_col).loc[sample_names]
    else:
        samples_metadata = pd.DataFrame(index=sample_names)
    
    if hasattr(args, 'dhs_annotations') and ('dist_tss' not in dhs_meta.columns):
        da = args.dhs_annotations
        if da is not None and os.path.exists(da):
            dhs_annotations = pd.read_table(da).set_index('dhs_id')
            dhs_meta['dist_tss'] = dhs_annotations['dist_tss']
    
    return parse_optional_args(args, mat, dhs_metadata=dhs_meta, samples_metadata=samples_metadata)


def mask_from_metadata(metadata: pd.DataFrame, column_name: str):
    if column_name is not None:
        return metadata[column_name].to_numpy().astype(bool)
    return np.ones(metadata.shape[0], dtype=bool)


def parse_optional_args(
        args,
        matrix: np.ndarray,
        samples_metadata: pd.DataFrame,
        dhs_metadata: pd.DataFrame
    ) -> NMFInputData:
    samples_m = read_mask(args.samples_mask)
    peaks_m = read_mask(args.peaks_mask)
    
    if args.samples_weights is not None or args.peaks_weights is not None:
        W_weights_vector = read_weights(
            args.samples_weights,
            matrix.shape[1],
            samples_metadata.index
        )#.astype(matrix.dtype)

        H_weights_vector = read_weights(
            args.peaks_weights,
            matrix.shape[0],
            dhs_metadata.index
        )#.astype(matrix.dtype)
    else:
        W_weights_vector = H_weights_vector = None
  
    return NMFInputData(
        matrix=matrix,
        samples_mask=samples_m,
        peaks_mask=peaks_m,
        samples_weights=W_weights_vector,
        peaks_weights=H_weights_vector,
        dhs_metadata=dhs_metadata,
        samples_metadata=samples_metadata,
    )

def read_mask(mask) -> np.ndarray:
    if mask is None:
        mask = np.ones(mask, dtype=bool)
    elif isinstance(mask, str):
        mask = np.loadtxt(mask, dtype=bool)
    return mask


def read_weights(weights_path, shape, sample_names=None) -> np.ndarray:
    if weights_path is not None:
        if os.path.splitext(weights_path)[-1] == 'npy':
            weights_vector = np.load(weights_path)
        else:
            weights_df = pd.read_table(weights_path).set_index('id')
            if sample_names is not None:
                weights_df = weights_df.loc[sample_names]
            weights_vector = weights_df['weight'].to_numpy().astype(dtype=dtype)
            assert np.all(np.isfinite(weights_vector)), 'Some provided weights are not finite!'
            assert weights_vector.shape[0] == shape
    else:
        weights_vector = np.ones(shape, dtype=dtype)
    
    return weights_vector / weights_vector.sum() * weights_vector.shape[0]


# Internal function to run NMF decomposition as part of pipeline
def main(nmf_input_data: NMFInputData, **extra_params):
    mat = nmf_input_data.matrix
    samples_m = nmf_input_data.samples_mask

    W_weights = nmf_input_data.samples_weights
    H_weights = nmf_input_data.peaks_weights


    samples_masked_matrix = mat[:, samples_m]
    print(samples_masked_matrix.shape, mat.shape)
    non_zero_rows = samples_masked_matrix.sum(axis=1) > 0
    print(non_zero_rows.sum())
    peaks_mask = nmf_input_data.peaks_mask & non_zero_rows

    matrix_samples_peaks_slice = samples_masked_matrix[peaks_mask, :]
    print(matrix_samples_peaks_slice.shape)

    if W_weights is not None or H_weights is not None:
        print('Using weighted NMF')
        W_weights_slice = W_weights[samples_m]
        H_weights_slice = H_weights[peaks_mask]
        is_weighted = True
    else:
        print('Using sklearn default NMF')
        W_weights_slice = H_weights_slice = None
        is_weighted = False

    model = initialize_model(
        n_components=args.n_components,
        extra_params=extra_params,
        is_weighted=is_weighted
    )
    
    print('Fitting NMF model', flush=True)
    W_np, H_np, model = run_NMF(
        model=model,
        X=matrix_samples_peaks_slice,
        W_weights=W_weights_slice,
        H_weights=H_weights_slice,
    )
    if args.project_masked_samples and samples_m.shape[0] > samples_m.sum():
        print('Projecting samples')
        if args.project_masked_peaks and (non_zero_rows.sum() != mat[:, samples_m].shape[0]):
            H_np = project_peaks(
                model,
                mat[:, samples_m][non_zero_rows, :],
                W_np,
                W_weights_slice,
                H_weights
            )
            W_np = project_samples(
                model,
                mat[non_zero_rows, :],
                H_np,
                W_weights,
                H_weights
            )
        else:
            W_np = project_samples(
                model,
                mat[peaks_mask, :],
                H_np,
                W_weights,
                H_weights_slice
            )

    return W_np, H_np, peaks_mask


def setup_parser():
    parser = argparse.ArgumentParser('Run NMF decomposition on a DHS x samples matrix')
    parser.add_argument('n_components', help='Number of components to use in NMF', type=int)
    parser.add_argument('prefix', help='Prefix for the output file')

    # Reading from npy matrix arguments
    parser.add_argument('--matrix', help='Path to matrix to run NMF on')
    parser.add_argument('--sample_names', help='Path to file with sample names')
    parser.add_argument('--dhs_meta', help='Path to DHS metadata file')
    parser.add_argument('--samples_mask', help='Mask of used samples, txt indicator', default=None)
    parser.add_argument('--peaks_mask', help='Mask of used peaks, txt indicator', default=None)

    # Reading from anndata argument
    parser.add_argument('--from_anndata', help='Path to AnnData file. If provided, ignore matrix, sample_names and dhs_meta fields.', default=None)
    parser.add_argument('--samples_mask_column', help='Column in samples metadata to use as mask', default=None)
    parser.add_argument('--dhs_mask_column', help='Column in dhs metadata to use as mask', default=None)

    # Optional arguments for both modes
    parser.add_argument('--project_masked_peaks', action='store_true', help='Project peaks for all masked samples', default=False)
    parser.add_argument('--project_masked_samples', action='store_true', help='Project peaks for all masked samples', default=False)
    parser.add_argument('--samples_weights', help='Path to samples weights (for weighted NMF)', default=None)
    parser.add_argument('--peaks_weights', help='Path to peaks weights (for weighted NMF)', default=None)
    parser.add_argument('--extra_params', help='Path to json with NMF params (overwrites default settings)', default=None)
    return parser


if __name__ == '__main__':
    parser = setup_parser()
    args = parser.parse_args()

    print('Parsing arguments')
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