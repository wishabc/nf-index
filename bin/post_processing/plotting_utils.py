import numpy as np
import pandas as pd
from tqdm import tqdm
import scipy.cluster.hierarchy as sch

def identity_agst(W, order=None):
    if order is None:
        order = np.arange(W.shape[0])
    elif order == 'agst':
        order = np.argsort((W / W.sum(axis=0)).mean(axis=1))
    order = np.array(order)
    assert len(order.shape) == 1 and order.shape[0] == W.shape[0]
    return np.tile(order, W.shape[1]).reshape(W.shape[::-1]).T


def sample_repr_dhs(valid_peaks_mask, n_to_sample=10_000):
    if valid_peaks_mask.sum() < n_to_sample:
        raise ValueError(f'Not enough valid peaks to sample: {valid_peaks_mask.sum()}')
    indices = np.arange(valid_peaks_mask.shape[0])[valid_peaks_mask]
    np.random.seed(0)
    return np.random.choice(indices, n_to_sample, replace=False) # indices


def get_sampled_dhs_indices(H, mask, n_to_sample=10_000, min_dhs_per_cluster=50, n_per_cluster=200, cluster_threshold=0.5):
    assert H.shape[1] == mask.shape[0]
    sub_mask = np.sort(sample_repr_dhs(mask, n_to_sample))
    H_subsampled = np.ascontiguousarray(H[:, sub_mask])
    
    order_s, clusters_s, linkage_s = hierarchical_cluster(H_subsampled, cluster_threshold=cluster_threshold)
    clusters, cluster_counts = np.unique(clusters_s, return_counts=True)
    selected_clusters = clusters[cluster_counts >= min_dhs_per_cluster]

    dhs_indices = []
    for s in selected_clusters:
        cluster_mask = clusters_s == s
        if cluster_mask.sum() < n_per_cluster:
            print(f'WARNING: not enough DHSs per cluster {cluster_mask.sum()}')
        sampled_from_mask_indices = sample_repr_dhs(cluster_mask, min(cluster_mask.sum(), n_per_cluster))
        full_dhs_index = sub_mask[order_s][np.sort(sampled_from_mask_indices)]
        # full_dhs_index = full_dhs_index[np.argsort(H_norm.max(axis=0)[full_dhs_index])]    
        dhs_indices.append(full_dhs_index)
    
    plotting_args = H_subsampled, order_s, clusters_s, linkage_s
    return plotting_args, dhs_indices


def hierarchical_cluster(W, cluster_threshold=0.7):
    distance_matrix = sch.distance.pdist((W).T, metric='cosine')
    linkage_matrix = sch.linkage(distance_matrix, method='average')
    dendrogram = sch.dendrogram(linkage_matrix, no_plot=True)
    leaves_order = np.array(dendrogram['leaves'])
    clusters = sch.fcluster(linkage_matrix, t=cluster_threshold, criterion='distance')
    
    return leaves_order, clusters[leaves_order], linkage_matrix


def construct_sample_labels(metadata):
    return metadata.apply(
            lambda row: f"{row['taxonomy_name']} {row['SPOT1_score']:.1f}",
            axis=1
        )


def average_profiles(matrix, metadata, columns, func=np.mean):
    gb = metadata.reset_index().groupby(columns)
    result = {}
    for g_name, g in tqdm(gb, total=len(gb.size())):
        result[g_name] = func(matrix[:, g.index], axis=1)

    return pd.DataFrame(result)


# Sets the first level annotation to hcl.
def set_column_0(columns, df):
    arrays = [columns]
    for i in range(1, df.columns.nlevels):
        arrays.append(df.columns.get_level_values(i))
    new_columns = pd.MultiIndex.from_arrays(arrays)
    return pd.DataFrame(df.values, columns=new_columns)