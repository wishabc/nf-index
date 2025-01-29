import numpy as np
import argparse
from tqdm import tqdm

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("suffix", help="Suffix for the file")
    parser.add_argument("samples_order", help="Samples order")
    parser.add_argument("outpath", help="Output file path")
    parser.add_argument("--dtype", default="float")

    args = parser.parse_args()

    dtype = dict(zip(('int', 'float', 'bool'), (np.int32, np.float32, bool)))[args.dtype]
    samples_order = np.loadtxt(args.samples_order, dtype=str)

    files = [f"{x}.{args.suffix}.txt" for x in samples_order]
    files = [np.loadtxt(file, dtype=dtype) for file in tqdm(files)]
    print('Concatenating vectors')
    matrix = np.stack(files)

    np.save(args.outpath, matrix)

