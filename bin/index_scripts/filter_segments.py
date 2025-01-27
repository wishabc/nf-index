import pandas as pd
import sys


assert len(sys.argv) == 2
df = pd.read_table(sys.argv[1]).query('success_fit & ~max_bg_reached')
df[['#chr', 'start', 'end']].to_csv(sys.stdout, sep='\t', index=False, header=None)