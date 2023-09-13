import pandas as pd
import sys

def main(df):
    df.groupby()
    
    
if __name__ == '__main__':
    df = pd.read_table(sys.argv[1], 
        names=['chr1', 'start1', 'end1', 'chunk_id', 'score', 
        'n_samples', 'n_peaks', 'n', 'summit', 'core_start', 'core_end',
        'chr2', 'start2', 'end2', 'AG', 'score', 'sum', 'summit']
    )
    main(df)