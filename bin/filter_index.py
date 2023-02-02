import pandas as pd
import numpy as np
import sys


index_cols = ["chr", "start", "end", "chunk_id", "score", "n_samples",
     "n_peaks", "dhs_width", "summit", "start_core", "end_core"]



def main(index, mask):
    index['avg_score'] = index.eval('score / n_samples')
    cutoff = index[mask]['avg_score'].quantile(0.05)
    print(f'Cutoff = {cutoff}')
    return (index['avg_score'] >= cutoff) * mask

if __name__ == '__main__':
    index = pd.read_table(sys.argv[1], header=None, names=index_cols)
    # mask - numpy array with zeros for peaks in ENCODE blacklisted regions
    mask = np.invert(np.loadtxt(sys.argv[2], dtype=bool))
    out_mask = main(index, mask)
    np.savetxt(sys.argv[3], out_mask)
