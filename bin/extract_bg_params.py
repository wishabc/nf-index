import argparse
import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Extract background parameters")
    parser.add_argument("bg_params_file", help="Path to the background params file overlapped with a bed file")
    parser.add_argument("bg_params_r", help="Path to the r parameters output file")
    parser.add_argument("bg_params_p", help="Path to the p parameters output file")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    data = pd.read_table(args.bg_params_file, na_values=".").drop_duplicates(
        subset=['#chr', 'start', 'end']
    )
    np.save(args.bg_params_r, data['bg_r'].values.astype(np.float32))
    np.save(args.bg_params_p, data['bg_p'].values.astype(np.float32))
