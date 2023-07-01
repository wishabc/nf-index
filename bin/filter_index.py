import pandas as pd
import numpy as np
import argparse


index_cols = ["#chr", "start", "end", "chunk_id", "score", "n_samples",
     "n_peaks", "dhs_width", "summit", "start_core", "end_core"]

autosomes = [f'chr{x}' for x in range(1, 23)]


def main(index, mask, include_lowsig_singletons=False):
    autosomes_mask = index['chr'].isin(autosomes)
    total_mask = autosomes_mask & mask
    if not include_lowsig_singletons:
        index['avg_score'] = index.eval('score / n_samples')
        cutoff = index[total_mask]['avg_score'].quantile(0.05)
        print(f'Cutoff = {cutoff}')
        return (index['avg_score'] >= cutoff) & total_mask
    return total_mask


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Filter index')
    parser.add_argument('index_file', help='Path to index_file')
    parser.add_argument('blacklisted_mask', help='Path to blacklisted_mask')
    parser.add_argument('output_mask', help='Path to directory to save output mask.')
    parser.add_argument('filtered_index', help='Path to directory to save filtered index.')
    parser.add_argument('--include_lowsig_singletons', help='Include low signal singletons',
        default=False, action="store_true")

    args = parser.parse_args()
    index_df = pd.read_table(args.index_file, header=None, names=index_cols)
    # mask - numpy array with zeros for peaks in ENCODE blacklisted regions
    mask = np.invert(np.loadtxt(args.blacklisted_mask, dtype=bool))
    out_mask = main(index_df, mask, args.include_lowsig_singletons)
    assert len(index_df.index) == len(out_mask)
    np.savetxt(args.output_mask, out_mask, fmt="%5i")
    index_df[index_cols][out_mask].to_csv(args.filtered_index, sep='\t', index=False)
