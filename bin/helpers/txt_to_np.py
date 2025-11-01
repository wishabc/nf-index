import numpy as np
import argparse



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input txt file")
    parser.add_argument("outpath", help="Output path to npy converted array")
    parser.add_argument("--dtype", default="float")

    args = parser.parse_args()

    dtype = dict(zip(('int', 'float', 'bool'), (np.int32, np.float32, bool)))[args.dtype]

    arr = np.loadtxt(args.input, dtype=dtype)
    np.save(args.outpath, arr)

