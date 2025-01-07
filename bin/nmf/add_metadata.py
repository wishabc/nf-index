import pandas as pd
import sys


def main(metadata, path_to_results):
    unique_prefix = metadata['prefix']
    metadata['W'] = path_to_results + '/' + unique_prefix + '/' + unique_prefix + '.W.npy'
    metadata['H'] = metadata['W'].str.replace('.W.npy', '.H.npy', regex=False)
    metadata['peaks_mask'] = metadata['W'].str.replace('.W.npy', '.non_zero_peaks_mask.txt', regex=False)
    return metadata


if __name__ == '__main__':
    metadata = pd.read_table(sys.argv[1]).drop_duplicates(subset='prefix')
    path_to_results = sys.argv[2]
    metadata = main(metadata, path_to_results)
    metadata.to_csv(sys.argv[3], sep='\t', index=False)
