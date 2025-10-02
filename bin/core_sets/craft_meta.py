import argparse
import pandas as pd
import numpy as np

def main():
    pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Get metadata for core sets"
    )
    parser.add_argument(
        "mapping",
        help="Path to the mapping file with term -> path_safe_id information"
    )
    parser.add_argument(
        "core_sets_outdir",
        help="Output directory where core sets are saved"
    )
    parser.add_argument(
        "output",
        help="Output file with core sets metadata"
    )
    parser.add_argument(
        "--fdrs",
        nargs="+",
        required=True,
        help="List of FDR thresholds"
    )
    args = parser.parse_args()

    mapping = pd.read_table(args.mapping)
    mapping['fdr'] = [args.fdrs] * len(mapping)
    mapping = mapping.explode('fdr')
    mapping['prefix'] = mapping['path_safe_id'].astype(str) + ".fdr" + mapping['fdr'].astype(str)
    mapping['base_path'] = args.core_sets_outdir + "/core_sets/" + mapping['fdr'].astype(str) + "/" + mapping['prefix']


    mapping['core_set_npy'] = mapping['base_path'] +  ".core_set.npy"
    mapping['core_set_bed'] = mapping['base_path'] + ".core_set.bed"
    mapping['saturation_curve'] = mapping['base_path'] +  ".saturation_curve.npy"
    mapping['saturation_curve_core'] = mapping['base_path'] + ".saturation_curve_core.npy"

    mapping['step_added'] = mapping['base_path'] + ".step_added.npy"
    mapping['mcv_by_step_stats'] = mapping['base_path'] + ".mcv_by_step_stats.npy"

    mapping['core_set_size'] = mapping['core_set_npy'].apply(
        lambda x: np.load(x).sum()
    )
    mapping = mapping.drop(columns=['base_path', 'prefix'])
    mapping.to_csv(args.output, sep='\t', index=False)


