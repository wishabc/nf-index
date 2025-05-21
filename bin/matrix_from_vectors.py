import numpy as np
import argparse
from tqdm import tqdm


def read_column(path, dtype):
    if path.endswith('.npy'):
        return np.load(path).astype(dtype)
    else:
        return np.loadtxt(path, dtype=dtype)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("suffix", help="Suffix for the file")
    parser.add_argument("samples_order", help="Samples order")
    parser.add_argument("outpath", help="Output file path")
    parser.add_argument("--dtype", default="float")
    parser.add_argument("--input_ext", default="txt", help="Input file extension")

    args = parser.parse_args()

    dtype = dict(zip(('int', 'float', 'bool'), (np.int32, np.float32, bool)))[args.dtype]
    samples_order = np.loadtxt(args.samples_order, dtype=str)

    files = [f"{x}.{args.suffix}.{args.input_ext}" for x in samples_order]
    
    files = [read_column(file, dtype=dtype) for file in tqdm(files)]
    print('Concatenating vectors')
    matrix = np.stack(files).T

    np.save(args.outpath, matrix)

