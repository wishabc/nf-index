import pandas as pd

import numpy as np

import sys



masterlist = pd.read_table(
    sys.argv[1],
    names=['chrom', 'start', 'end', 'dhs_id', 'summit']
)

sample_id = sys.argv[2]


jitters = np.load(sys.argv[3], mmap_mode='r')

jitter_offsets = jitters[sample_id]

masterlist['summit'] += jitter_offsets.astype(int)

masterlist['start'] = masterlist['summit'].astype(int)
masterlist['end'] = masterlist['start'] + 1

masterlist = masterlist.sort_values(['chrom', 'start', 'end'])

masterlist[['chrom', 'start', 'end', 'dhs_id', 'summit']].to_csv(
    sys.argv[4],
    sep='\t',
    header=False,
    index=False,
)

np.savetxt(sys.argv[5], masterlist.index.to_numpy(), fmt="%d")