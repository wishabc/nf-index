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


def add_normalization_params(adata, params):
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
            loaded_params = np.load(param)
            for key in loaded_params:
                adata.uns[f'lowess_normalization_params.{key}'] = loaded_params[key]
        else:
            raise ValueError(f'Unknown parameter file type: {param}')


def main(
        adata,
        matrices,
        params,
        variance_partition_results,
        mask
    ):
    matrices_types = [
        'scale_factors.mean_normalized',
        'vst'
    ]
    matrices_mapping = get_matrices_mapping_by_types(matrices, matrices_types)

    add_matrices_to_anndata(adata, matrices_mapping, mask)
    add_normalization_params(adata, params)
    
    adata.varm['variance_partition'] = expand_vp_results(
        variance_partition_results,
        adata.var.index,
        mask
    ).values
    return adata


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Merge results of normalization into existing AnnData object."
    )
    parser.add_argument("input", help="Input Zarr file")
    parser.add_argument("output", help="Output Zarr file")
    parser.add_argument("variance_partition_result", help="Variance partition annotated masterlist bed")
    parser.add_argument("deseq_design_formula", help="Formula used for variance partition")
    parser.add_argument("variance_partition_formula", help="Formula used for variance partition")
    
    parser.add_argument("--normalization_layer", help="Name of the layer used for normalization", default="counts")
    parser.add_argument("--dhs_mask_name", help="DHS mask column name")
    parser.add_argument(
        "--matrices",
        nargs="+",
        required=True,
        help="Matrices to add to the AnnData object"
    )
    parser.add_argument(
        "--params",
        nargs="+",
        default=[],
        help="Normalization parameters to reproduce normalization (RDS, json, npz files)"
    )
    args = parser.parse_args()

    adata = read_zarr_backed(args.input)

    variance_partition_results = pd.read_table(args.variance_partition_result)
    variance_partition_results = variance_partition_results[
        [x for x in variance_partition_results.columns if not x.startswith('V')]
    ]

    mask = get_mask_from_column_name(adata, args.dhs_mask_name)

    adata = main(
        adata,
        args.matrices,
        args.params,
        variance_partition_results,
        mask=mask
    )
    adata.uns['normalization_layer'] = args.normalization_layer
    adata.uns['variance_partition_formula'] = args.variance_partition_formula
    adata.uns['deseq_design_formula'] = args.deseq_design_formula

    adata.write_zarr(args.output)