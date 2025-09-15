import numpy as np
import pandas as pd
import argparse
import json
from genome_tools.data.anndata import read_zarr_backed
import os

import scipy.sparse as sp
import dataclasses
from nmf_tools.nmf import NMFModel, read_weights

dtype = np.float64


@dataclasses.dataclass
class NMFInputData:
    matrix: sp.csr_matrix
    samples_mask: np.ndarray
    peaks_mask: np.ndarray
    samples_weights: np.ndarray
    peaks_weights: np.ndarray
    samples_metadata: pd.DataFrame
    dhs_metadata: pd.DataFrame
    mode: str



##################
### Parse args ###
##################

def parse_nmf_args(args) -> NMFInputData:
    if args.from_anndata is not None:
        return parse_args_anndata(args)
    return parse_args_matrix(args)


def parse_args_anndata(args):
    print('Reading AnnData')
    adata = read_zarr_backed(args.from_anndata)

    # ### Debug part start ###
    # indices = np.arange(adata.shape[1])
    # np.random.seed(42)
    # selected_indices = indices[np.random.rand(adata.shape[1]) < 0.05]
    # print(f'adata shape: {adata.shape[1]}')
    # selected_indices = np.sort(selected_indices)
    # np.save('/home/sabramov/tmp/selected_indices.npy', selected_indices)
    # adata = adata[:, selected_indices]
    # ### Debug part end ###

    matrix = adata.layers['binary'].T

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
    
    return parse_optional_args(
        args, mat,
        dhs_metadata=dhs_meta,
        samples_metadata=samples_metadata
    )


def mask_from_metadata(metadata: pd.DataFrame, column_name: str):
    if column_name is not None:
        return metadata[column_name].to_numpy().astype(bool)
    return np.ones(metadata.shape[0], dtype=bool)


def parse_optional_args(
        args,
        matrix: sp.csr_matrix,
        samples_metadata: pd.DataFrame,
        dhs_metadata: pd.DataFrame
    ) -> NMFInputData:
    samples_m = read_mask(args.samples_mask)
    peaks_m = read_mask(args.peaks_mask)
    
    W_weights_vector = H_weights_vector = None
    if args.samples_weights is not None:
        W_weights_vector = read_weights(
            args.samples_weights,
            names=samples_metadata.index
        ).values


    if args.peaks_weights is not None:
        H_weights_vector = read_weights(
            args.peaks_weights,
            names=dhs_metadata.index
        ).values

    return NMFInputData(
        matrix=matrix,
        samples_mask=samples_m,
        peaks_mask=peaks_m,
        samples_weights=W_weights_vector,
        peaks_weights=H_weights_vector,
        dhs_metadata=dhs_metadata,
        samples_metadata=samples_metadata,
        mode=args.mode
    )

def read_mask(mask) -> np.ndarray:
    if mask is None:
        mask = np.ones(mask, dtype=bool)
    elif isinstance(mask, str):
        mask = np.loadtxt(mask, dtype=bool)
    return mask



# Internal function to run NMF decomposition as part of pipeline
def main(nmf_input_data: NMFInputData, n_components: int, **extra_params):
    mat = nmf_input_data.matrix
    samples_mask = nmf_input_data.samples_mask
    W_weights = nmf_input_data.samples_weights
    H_weights = nmf_input_data.peaks_weights

    # Slice matrix by masked samples
    samples_matrix = mat[:, samples_mask]
    print(samples_matrix.shape, mat.shape)

    # Keep peaks with non-zero signal in selected samples
    non_zero_rows = samples_matrix.sum(axis=1) > 0
    print(non_zero_rows.sum())
    peaks_mask = nmf_input_data.peaks_mask & non_zero_rows

    matrix_slice = samples_matrix[peaks_mask, :]

    # Handle weights
    if W_weights is not None or H_weights is not None:
        W_weights_slice = W_weights[samples_mask] if W_weights is not None else None
        H_weights_slice = H_weights[peaks_mask] if H_weights is not None else None
    else:
        print("Using sklearn default NMF")
        W_weights_slice = H_weights_slice = None

    # Initialize NMF model
    nmf_model = NMFModel(n_components=n_components, extra_params=extra_params)

    print("Fitting NMF model", flush=True)
    if nmf_input_data.mode == "weighted":
        W_np, H_np = nmf_model.fit_transform(
            X=matrix_slice,
            W_weights=W_weights_slice,
            H_weights=H_weights_slice,
        )
    else:
        X = matrix_slice.copy()
        if W_weights_slice is not None:
            X = X.multiply(W_weights_slice[:, None] ** 0.5)
        if H_weights_slice is not None:
            X = X.multiply(H_weights_slice[None, :] ** 0.5)

        W_np, H_np = nmf_model.fit_transform(X=X)

    # Optionally project back masked samples
    if args.project_masked_samples and samples_mask.shape[0] > samples_mask.sum():
        print("Projecting samples")
        W_np = nmf_model.project_samples(
            mat[peaks_mask, :],
            H_np,
            W_weights=W_weights,
            H_weights=H_weights_slice,
        )
        W_weights_slice = W_weights

    #  adjust output for modified_X mode
    if nmf_input_data.mode == "modified_X":
        return (
            W_np / np.sqrt(W_weights_slice[None, :]),
            H_np / np.sqrt(H_weights_slice[:, None]),
            peaks_mask,
        )

    return W_np, H_np, peaks_mask



def setup_parser():
    parser = argparse.ArgumentParser('Run NMF decomposition on a DHS x samples matrix')
    parser.add_argument('n_components', help='Number of components to use in NMF', type=int)
    parser.add_argument('prefix', help='Prefix for the output file')

    # Reading from npy matrix arguments, deprecated
    parser.add_argument('--matrix', help='Path to matrix to run NMF on')
    parser.add_argument('--sample_names', help='Path to file with sample names')
    parser.add_argument('--dhs_meta', help='Path to DHS metadata file')
    parser.add_argument('--samples_mask', help='Mask of used samples, txt indicator', default=None)
    parser.add_argument('--peaks_mask', help='Mask of used peaks, txt indicator', default=None)

    # Reading from anndata argument
    parser.add_argument('--from_anndata', help='Path to AnnData file. If provided, ignore matrix, sample_names and dhs_meta fields.', default=None)
    parser.add_argument('--samples_mask_eval', help='Eval expression on the samples metadata to use as samples mask', default=None)
    parser.add_argument('--dhs_mask_eval', help='Eval expression on the dhs metadata to use as DHS mask', default=None)

    # Optional arguments for both modes
    parser.add_argument('--mode', choices=('weighted', 'modified_X'), help='Whether use weighted NMF or update X as sqrt(weight) * X and then apply inverse transform for H and W', default=None)
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

    W_np, H_np, peaks_mask = main(nmf_data, args.n_components, **extra_params)

    print('Saving results')
    np.save(f'{args.prefix}.W.npy', W_np) # samples x components
    np.save(f'{args.prefix}.H.npy', H_np.T) # peaks x components
    np.savetxt(f'{args.prefix}.non_zero_peaks_mask.txt', peaks_mask, fmt="%d")
    np.savetxt(f'{args.prefix}.samples_mask.txt', nmf_data.samples_mask, fmt="%d")
    #np.savetxt(f'{args.prefix}.weights_vector.txt', weights_vector)