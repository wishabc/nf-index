#!/bin/usr/python

import sys
import logging
import datatable as dt
import argparse
import numpy as np
import scipy.sparse as sp


logger = logging.getLogger(__name__)
handler = logging.StreamHandler(sys.stdout)
logger.setLevel('INFO')
formatter = logging.Formatter('%(asctime)s  %(levelname)s  %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)


def read_matrix(input_file, dtype=int):
    # Make sure to say there is NO header in this file. Otherwise will be off by 1
    df = dt.fread(input_file, header=False)
    logger.info('Converting to numpy array')
    return df.to_numpy().astype(dtype)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Converts matrix in txt format to csr')
    parser.add_argument('matrix', help='Path to matrix file')
    parser.add_argument('outpath', help='Path to output binary matrix file with .npy extension')
    parser.add_argument('--dtype', choices=('int', 'float', 'bool'), help='dtype of the matrix file', default='float')
    parser.add_argument('--mask', help='Mask file to filter the matrix by rows', default=None)
    parser.add_argument('--fmt', help='Format to save matrix in', choices=('dense', 'sparse'), default='dense')
    args = parser.parse_args()
    input_path = args.matrix
    out_path = args.outpath
    logger.info('Starting processing')

    dtype = dict(zip(('int', 'float', 'bool'), (np.uint32, np.float32, bool)))[args.dtype]
    matrix_dense = read_matrix(input_path, dtype=dtype)
    logger.info(f'Matrix size: {matrix_dense.shape}. '
                f'Density: {(matrix_dense != 0).sum() / matrix_dense.size}'
                )
    if args.mask is not None:
        mask = np.loadtxt(args.mask, dtype=bool)
        assert mask.shape[0] == matrix_dense.shape[0]
        matrix_dense = matrix_dense[mask, :]
    if args.fmt == 'dense':
        np.save(args.outpath, matrix_dense)
    elif args.fmt == 'sparse':
        matrix_sparse = sp.coo_matrix(matrix_dense)
        sp.save_npz(out_path, matrix_sparse)
