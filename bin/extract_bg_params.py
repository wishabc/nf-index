import argparse
import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Extract background parameters")
    parser.add_argument("bg_params_file", help="Path to the bg params file overlapped with a bed file")
    parser.add_argument("output", help="Path to the output file")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    data = pd.read_table(args.bg_params_file, na_values=".")
    np.save(args.output, data[['bg_r', 'bg_p']].values.astype(np.float32))
