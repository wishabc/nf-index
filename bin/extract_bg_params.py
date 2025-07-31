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

    header = np.loadtxt(args.header, dtype=str)
    data = pd.read_table(sys.stdin, header=None, names=['#chr1', 'start1', 'end1', *header])
    np.save(args.output, data[['bg_r', 'bg_p']].values.astype(np.float32))