import numpy as np
import sys
from tqdm import tqdm

if __name__ == '__main__':
    # TODO: wrap in a function
    chunk_size = int(sys.argv[1])
    matrix = np.load(sys.argv[2])
    print(matrix.shape)
    for start_index in tqdm(list(range(0, matrix.shape[0], chunk_size))):
        end_index = start_index + chunk_size
        np.savetxt(
            f"{sys.argv[3]}.{start_index}",
            matrix[start_index:end_index, :]
        )
