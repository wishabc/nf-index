import sys
import pandas as pd
import numpy as np
import os
from helpers import add_matrices_to_anndata, get_matrices_mapping_by_types
from genome_tools.data.anndata import read_zarr_backed
import json
import base64
import argparse

from helpers import get_mask_from_column_name


def expand_vp_results(variance_partition_results: pd.DataFrame, index, mask):
    expanded_variance_partition_results = pd.DataFrame(
        data=np.nan,
        index=index,
        columns=variance_partition_results.columns
    )
    expanded_variance_partition_results.loc[mask, :] = variance_partition_results.values
    return expanded_variance_partition_results


def add_normalization_params(adata, params, mask=None):
    for param in params:
        if param.endswith('.dispersion_function.RDS'):
            with open(param, 'rb') as f:
                adata.uns['vst_dispersion_function_rds.b64_encoded'] = base64.b64encode(
                    f.read()
                ).decode('utf-8')
        elif param.endswith('json'):
            with open(param) as f:
                adata.uns['lowess_normalization_params'] = json.load(f)
        elif param.endswith('npz'):
            if mask is None:
                mask = np.ones(adata.n_vars, dtype=bool)
            loaded_params = np.load(param)
            for key, array in loaded_params.items():
                if key == "norm_factors_geometric_mean":
                    result = np.full(adata.n_vars, np.nan, dtype=array.dtype)
                    result[mask] = array
                    adata.varm[f'lowess_normalization_params.{key}'] = result
                else:
                    adata.uns[f'lowess_normalization_params.{key}'] = array
        else:
            raise ValueError(f'Unknown parameter file type: {param}')


def add_bed_files_to_anndata(adata, bed_files, mask):

    for bed_file in bed_files:
        bed_df = pd.read_table(bed_file)
        if 'vp_annotated' in bed_file:
            key = 'variance_partition'
            bed_df = bed_df[[x for x in bed_df.columns if not x.startswith('V')]]
        else:
            raise ValueError(f'Unknown bed file type: {bed_file}')

        adata.varm[key] = expand_vp_results(
            bed_df,
            adata.var.index,
            mask
        )

def main(
        adata,
        matrices,
        params,
        bed,
        mask
    ):
    expected_matrices_type = [
        'scale_factors.mean_normalized',
        'vst',
        'bg_corrected_agg_cutcounts',
    ] + list(adata.layers.keys())
    matrices_mapping = get_matrices_mapping_by_types(matrices, expected_matrices_type)

    add_matrices_to_anndata(adata, matrices_mapping, mask)
    add_normalization_params(adata, params, mask)
    
    add_bed_files_to_anndata(adata, bed, mask)
    return adata


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Merge results of normalization into existing AnnData object."
    )
    parser.add_argument("input", help="Input Zarr file")
    parser.add_argument("output", help="Output Zarr file")
    parser.add_argument("--dhs_mask_name", help="DHS mask column name")
    parser.add_argument(
        "--layers",
        nargs="+",
        required=True,
        help="Matrices to add to the AnnData layers"
    )
    parser.add_argument(
        "--varm",
        help="BED files to add to the AnnData varm",
        nargs="*",
        default=[],
    )
    parser.add_argument(
        "--params",
        nargs="+",
        default=[],
        help="Normalization parameters to reproduce normalization (RDS, json, npz files)"
    )
    parser.add_argument(
        "--uns",
        nargs="+",
        default=[],
        help="Strings to save in uns (e.g. R formulae). key=value format"
    )
    args = parser.parse_args()

    adata = read_zarr_backed(args.input)

    mask = get_mask_from_column_name(adata, args.dhs_mask_name)
    adata.uns['dhs_mask_name'] = args.dhs_mask_name

    adata = main(
        adata,
        args.layers,
        args.params,
        args.varm,
        mask=mask
    )
    for arg in args.uns:
        key, value = arg.split('=', 1)
        adata.uns[key] = value

    adata.write_zarr(args.output)