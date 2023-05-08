import os

def get_matrix_path(dir_path_mat, method, data_slice_size, n_components, prefix):
    return os.path.join(dir_path_mat, f'{method}.{data_slice_size}.{n_components}.{prefix}_matrix.npy')