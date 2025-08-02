"""
       | python3 $moduleDir/bin/extract_bg_params.py 
            --header header.txt 
            --output ${name}
"""

import argparse
import numpy as np
import pandas as pd
import sys


def parse_args():
    parser = argparse.ArgumentParser(description="Extract background parameters")
    parser.add_argument("--header", required=True, help="Path to the header file")
    parser.add_argument("--output", required=True, help="Path to the output file")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    header = pd.read_table(args.header).columns
    print(header)
    names = ['#chr1', 'start1', 'end1', 'n_overlap', *header]
    print(names)
    data = pd.read_table(sys.stdin, header=None, names=names, na_values=".")
    print(data.head())
    np.save(args.output, data[['bg_r', 'bg_p']].values.astype(np.float32))