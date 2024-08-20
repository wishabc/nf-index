import pandas as pd
import sys


metadata = pd.read_table(sys.argv[1])
path_to_results = sys.argv[2]
metadata['W'] = path_to_results + '/' + metadata['prefix'] + '/' + metadata['prefix'] + '.' + metadata['n_components'].astype(str) + '.W.npy'
metadata['H'] = metadata['W'].str.replace('.W.npy', '.H.npy', regex=False)

metadata.to_csv(sys.argv[3], sep='\t', index=False)