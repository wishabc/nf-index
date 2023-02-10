import pandas as pd
import numpy as np
import sys


index_cols = ["chr", "start", "end", "chunk_id", "score", "n_samples",
     "n_peaks", "dhs_width", "summit", "start_core", "end_core"]

autosomes = [f'chr{x}' for x in range(1, 23)]


def main(index, mask):
    autosomes_mask = index['chr'].isin(autosomes)
    total_mask = autosomes_mask & mask
    index['avg_score'] = index.eval('score / n_samples')
    cutoff = index[total_mask]['avg_score'].quantile(0.05)
    print(f'Cutoff = {cutoff}')
    return (index['avg_score'] >= cutoff) & total_mask


if __name__ == '__main__':
    index_df = pd.read_table(sys.argv[1], header=None, names=index_cols)
    # mask - numpy array with zeros for peaks in ENCODE blacklisted regions
    mask = np.invert(np.loadtxt(sys.argv[2], dtype=bool))
    out_mask = main(index_df, mask)
    assert len(index_df.index) == len(out_mask)
    np.savetxt(sys.argv[3], out_mask, fmt="%5i")
    index_df[out_mask].to_csv(sys.argv[4], sep='\t', index=False, header=None)
