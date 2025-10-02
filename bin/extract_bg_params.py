import argparse
import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Extract background parameters")
    parser.add_argument("bg_params_file", help="Path to the background params file overlapped with a bed file")
    parser.add_argument("masterlist_file", help="Path to the masterlist bed file")
    parser.add_argument("mean_background", help="Path to the r parameters output file")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    data = pd.read_table(args.bg_params_file, na_values=".")
    dhs_ids = pd.read_table(args.masterlist_file, names=['#chr', 'start', 'end', 'dhs_id', 'summit'])['dhs_id'].values
    data = data.set_index('dhs_id').loc[dhs_ids, :]
    data['mean_bg'] = data['bg_r'] * data['bg_p'] / (1 - data['bg_p'])
    np.save(args.mean_background, data['mean_bg'].values.astype(np.float32))
