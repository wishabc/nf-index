import os

def get_matrix_path(dir_path_mat, method, data_slice_size, n_components, mat):
    return os.path.join(dir_path_mat, f'{method}.{data_slice_size}.{n_components}.{mat}_matrix.npy')