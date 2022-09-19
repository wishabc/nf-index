import sys
import datatable as dt
import logging
import argparse
import numpy as np

logger = logging.getLogger(__name__)
handler = logging.StreamHandler(sys.stdout)
logger.setLevel('INFO')
formatter = logging.Formatter('%(asctime)s  %(levelname)s  %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)


def main(raw_tags_path, normed_path, scale_factor_path):
    logger.info('Reading raw tags...')
    raw_tags = np.load(raw_tags_path)

    logger.info('Reading new normed matrix ...')
    as_normed_matrix = np.load(normed_path)

    logger.info('Calculating SF ...')
    sf = raw_tags / as_normed_matrix
    i = (raw_tags == 0) | (as_normed_matrix == 0)
    sf[i] = 1

    sf[np.isnan(sf)] = 1
    sf[~np.isfinite(sf)] = 1

    sf_geomean = np.exp(np.mean(np.log(sf), axis=1))
    as_scale_factors = (sf.T / sf_geomean).T
    np.save(scale_factor_path, as_scale_factors)

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Find scale factors of the matrix')
    parser.add_argument('signal_matrix', help='Path to matrix with read counts for each peak in every sample')
    parser.add_argument('normed_matrix', help='Path to normed matrix')
    parser.add_argument('output', help='File to save scaled factors into.')
    p_args = parser.parse_args()
    main(raw_tags_path=p_args.signal_matrix, 
         normed_path=p_args.normed_matrix,
         scale_factor_path=p_args.output)