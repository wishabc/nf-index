import h5py
import numpy as np
import argparse


def main(binary_matrix, vst_matrix, sample_names, out_path):
    samples = samples = np.char.encode(sample_names, 'ascii').astype('S')
    with h5py.File(out_path,  'w') as hf:
        hf.create_dataset('vst', data=vst_matrix.astype(np.float32))
        if binary_matrix is not None:
            hf.create_dataset('binary', data=binary_matrix.astype(bool))
        hf.create_dataset('sample_names', data=samples)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Store matrices in h5 format')
    parser.add_argument('vst', help='Path to vst normalized matrix')
    parser.add_argument('sample_names', help='Path to file with names of the columns of vst/binary matrix')
    parser.add_argument('outpath', help='Path to output h5 file')
    parser.add_argument('--binary', help='Path to binary matrix. Rows and columns are expected to be in the same order as in vst matrix', default=None)
    args = parser.parse_args()

    binary = args['binary'] if args['binary'] is None else np.load(args['binary'])
    vst = np.load(args['vst'])
    sample_names = np.loadtxt(args['sample_names'], dtype=str)
    main(binary, vst, sample_names, args['outpath'])