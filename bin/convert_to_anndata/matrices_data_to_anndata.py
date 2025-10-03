import sys
import pandas as pd
from helpers import add_matrices_to_anndata, get_matrices_mapping_by_types
from genome_tools.data.anndata import read_zarr_backed
import anndata as ad


def main(adata, matrices):
    matrices_types = [
        'counts', 'binary', 'density',
        'mean_bg_agg_cutcounts', 'neglog10_pvals'
    ] + list(adata.layers.keys())
    matrices_mapping = get_matrices_mapping_by_types(matrices, matrices_types)

    add_matrices_to_anndata(adata, matrices_mapping)

    adata.obs['n_peaks'] = adata.layers['binary'].sum(axis=1).A1
    return adata


if __name__ == '__main__':
    adata_obj = read_zarr_backed(sys.argv[1])
    samples_meta = pd.read_table(sys.argv[2], low_memory=False).set_index('sample_id')
    samples_meta[
        samples_meta.select_dtypes(include=['object', 'O', 'category']).columns
    ] = samples_meta.select_dtypes(include=['object', 'O', 'category']).fillna('None')

    
    missing_indices = adata_obj.obs.index.difference(samples_meta.index)
    if len(missing_indices) > 0:
        print(f"Warning: Missing indices: {missing_indices}")
        print("Creating a new anndata object (dropping X and obs)")
        adata_obj = ad.AnnData(X=None, obs=samples_meta, var=adata_obj.var)
    else:
        samples_meta = samples_meta.loc[adata_obj.obs.index, :]
        overlap = samples_meta.columns.intersection(adata_obj.obs.columns)

        if len(overlap) > 0:
            for col in overlap:
                assert samples_meta[col].equals(adata_obj.obs[col]), f"Column '{col}' differs between samples_meta and adata.obs"

        adata_obj.obs = samples_meta.drop(columns=overlap).join(adata_obj.obs)


    matrices = sys.argv[4:]

    main(adata_obj, matrices).write_zarr(sys.argv[3])
