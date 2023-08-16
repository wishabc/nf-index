import numpy as np
import sys


if __name__ == '__main__':
    # TODO: wrap in function
    chunk_size = int(sys.argv[1])
    matrix = np.load(sys.argv[2])
    for start_index in range(0, matrix.shape[0], chunk_size):
        end_index = start_index + chunk_size
        np.savetxt(
            f"{sys.argv[3]}.{start_index}",
            matrix[start_index:end_index, :]
        )
