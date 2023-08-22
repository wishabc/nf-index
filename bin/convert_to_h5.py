import h5py
import numpy as np
import sys


def main(binary_matrix, vst_matrix, sample_names, out_path):
    samples = samples = np.char.encode(sample_names, 'ascii').astype('S')
    with h5py.File(out_path,  'w') as hf:
        hf.create_dataset('vst', data=vst_matrix.astype(np.float32))
        hf.create_dataset('binary', data=binary_matrix.astype(bool))
        hf.create_dataset('sample_names', data=samples)


if __name__ == '__main__':
    binary = np.load(sys.argv[1])
    vst = np.load(sys.argv[2])
    sample_names = np.loadtxt(sys.argv[3], dtype=str)
    main(binary, vst, sample_names, sys.argv[4])