import pyBigWig
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import numpy as np
import pandas as pd
import sys

def extract_data(bw_path, interval):
    with pyBigWig.open(bw_path) as bw:
        data = bw.values(interval[0], interval[1], interval[2], numpy=True)
    return data

def main(bw_paths, interval, funcs):

    # Initialize the result array based on the number of BigWig files
    result = np.empty((len(bw_paths), interval[2] - interval[1]))

    # Process the BigWig files in parallel
    with ThreadPoolExecutor(max_workers=20) as executor:
        futures = {
            executor.submit(extract_data, bw_path, interval): i 
            for i, bw_path in enumerate(bw_paths)
        }
        for future in tqdm(as_completed(futures), total=len(futures)):
            i = futures[future]
            result[i, :] = future.result()
    
    results = {fname: func(result, axis=0) for fname, func in funcs.items()}
    return results


def write_to_bed(filename, interval, values):
    with open(filename, 'w') as bed_file:
        for i, value in enumerate(values):
            bed_file.write(f"{interval[0]}\t{interval[1] + i}\t{interval[1] + i + 1}\t{value}\n")


if __name__ == '__main__':
    met = pd.read_table(sys.argv[2])
    bw_paths = met['normalized_density_bw'].values
    str_interval = sys.argv[1]
    chrom, start, end = str_interval.split('_')
    interval = (chrom, int(start), int(end))
    funcs = {'mean': np.mean, 'std': np.std, 'median': np.median, 'min': np.min, 'max': np.max}
    results = main(bw_paths, interval, funcs)
    for fname, data in results.items():
        write_to_bed(f'{fname}.{str_interval}.bed', interval, data)