import sys
import numpy as np
import matplotlib.pyplot as plt
import anndata as ad
import scanpy as sc
import seaborn as sns
from scipy.sparse import csr_matrix
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist, squareform
from numba import jit


def pairwise_distances(X, metric='euclidean'):
    return squareform(pdist(X.T, metric))


def get_entropy_same_num(samples_meta):
    same = []
    entropy_same = dict()
    for s in same:
        for x in s:
            if x in entropy_same:
                entropy_same[x] += list(s - {x})
            else:
                entropy_same[x] = list(s - {x})

    validation_idx = samples_meta['core_annotation2'].notna()

    labels_alphabetical = list(np.unique(samples_meta.loc[validation_idx, 'core_annotation2'].to_numpy()))
    return {labels_alphabetical.index(k): [labels_alphabetical.index(x) for x in v] for k, v in entropy_same.items()}


def calc_entropy(euclid_dist, samples_meta, entropy_same_num):
    validation_idx = samples_meta['core_annotation2'].notna()

    entropy = np.full(validation_idx.shape, np.nan)
    entropy[validation_idx] = get_entropy_scores(
        euclid_dist[validation_idx, :][:, validation_idx], 
        np.unique(samples_meta.loc[validation_idx, 'core_annotation2'], return_inverse=True)[1],
        entropy_same_num
    )

    return entropy

def get_entropy_scores(euclid, labels_array, same_num_dict=None):
    annotation_entropy_scores = []
    for i in range(len(labels_array)):
        annotation_entropy_scores.append(get_metric(euclid[i, :], labels_array, same_num_dict=same_num_dict))
    
    annotation_entropy_scores = np.array(annotation_entropy_scores)
    
    return annotation_entropy_scores

def get_metric(distances, labels, same_num_dict=None):
    n = len(distances)
    assert len(labels) == n
    
    sorted_indices = np.argsort(distances)
    sorted_labels = labels[sorted_indices]

    n_nearest = (sorted_labels == sorted_labels[0]).sum() * 0.5
    p = np.power(0.1, 1/n_nearest)
    
    p_values = np.power(p, np.arange(n))
    p_values /= p_values.sum()

    # Ignore similar annotations in entropy score, consider them same as original
    if same_num_dict is not None and sorted_labels[0] in same_num_dict.keys():
        for v in same_num_dict[sorted_labels[0]]:
            sorted_labels[sorted_labels == v] = sorted_labels[0]
    
    
    label_count = accumulate_counts(sorted_labels, p_values, n)
    
    values = np.array(list(label_count.values()))
    values = values[values != 0]
    entropy = -np.sum(values * np.log2(values))
    
    return entropy

@jit(nopython=True)
def accumulate_counts(sorted_labels, p_values, n):
    label_count = {}
    for i in range(n):
        label = sorted_labels[i]
        if label not in label_count:
            label_count[label] = 0.0
        label_count[label] += p_values[i]
    return label_count


class ModelPlotter():
    def __init__(self, adata, nearest_neighbors, resolution,
                 na_color='#707070',  comp_prop=0.1, random_state=234234):
        self.nearest_neighbors = nearest_neighbors
        self.resolution = resolution
        self.comp_prop = comp_prop
        self.random_state = random_state

        self.adata = adata

        self.base_columns = ['core_system', 'super_state', 'growth_stage', 'raj_dataset', 'core_annotation2', 'extended_annotation2',
                            'treatment_class']
        self.other_columns = ['subsystem', 'donor_sex', 'disease_pathology', 'sampling_origin', 'organ', 'organ_region',
            'cell_type',
            'cell_line_type',
            'tissue_type',
             'raj_dataset', 'state', 'tag',
            'core_ontology_term',
            'frac_method',
            'num_peaks_downsampled',
            'ln_finished_date'
            ]

        self.data = self.adata.obsm['embedding'].T

        if  self.adata.obs['ln_finished_date'].dtype == 'object':
            self.adata.obs['ln_finished_date'] = self.adata.obs['ln_finished_date'].str.slice(0, 7).apply(self.get_date_num)

        self.base_columns += ['major_component']
        self.adata.obs['major_component'] = self.get_labels_from_nmf(self.data)

        self.na_color = na_color
        self.colormaps = {}

    def get_labels_from_nmf(self, data):
        argsort = np.argsort(np.array(data), axis=0)[::-1]
        cumul_sorted = np.cumsum(np.sort(np.array(data) / np.array(data).sum(axis=0), axis=0)[::-1], axis=0)
        idxs = np.concatenate([np.ones((1, np.array(data).shape[1]), dtype=bool), cumul_sorted < self.comp_prop], axis=0)[:-1]
        return np.array([','.join(map(str, np.sort(a[b]))) for a, b in zip(argsort.T, idxs.T)])


    def set_colormaps(self):
        for column in self.base_columns + self.other_columns + ['leiden', 'louvain']:
            mapped, cmp = get_discrete_colormap(self.adata.obs[column], na_color=self.na_color)
            self.colormaps[column] = {
                'colormap': cmp,
                'mapped_colors': mapped,
            }

    def set_distances(self):
        #Init adata fields with something
        # sc.pp.neighbors(self.adata, n_neighbors=5, use_rep='X', random_state=self.random_state)

        #Then fill them with distances (not sure if need to truncate to 0 at nn)
        self.adata.obsp['connectivities'] = csr_matrix(1 - self.adata.obsp['distance_matrix'])
        self.adata.obsp['distances'] = csr_matrix(self.adata.obsp['distance_matrix'])

    def calculate_clustering(self):
        sc.pp.neighbors(self.adata, n_neighbors=self.nearest_neighbors, use_rep='X', random_state=234234)
        sc.tl.leiden(self.adata, resolution=self.resolution)
        sc.tl.louvain(self.adata, resolution=self.resolution)

    def calculate_tsne(self):
        sc.tl.tsne(self.adata, use_rep='embedding')
        self.adata.obs.loc[:, ['TSNE1', 'TSNE2']] = self.adata.obsm['X_tsne']

    def average_distance_matrix(self):
        work_data = self.adata.obs[self.adata.obs['core_annotation2'].notna()]
        distances = self.adata.obsp['distance_matrix'][:, self.adata.obs['core_annotation2'].notna()][self.adata.obs['core_annotation2'].notna(), :]
        unique_labels = np.unique(work_data['core_annotation2'])
        M = len(unique_labels)
        
        groupped_distance_matrix = np.zeros((M, M), dtype=np.float_)
        for i in range(M):
            idx_i = work_data['core_annotation2'] == unique_labels[i]
            matrix_i = distances[idx_i, :]
            n = idx_i.sum()
            for j in range(i, M):
                idx_j =  work_data['core_annotation2'] == unique_labels[j]
                matrix_ij = matrix_i[:, idx_j]
                
                if i == j:
                    groupped_distance_matrix[i, j] = (np.sum(matrix_ij) - np.trace(matrix_ij)) / (n * (n - 1))
                else:
                    groupped_distance_matrix[i, j] = np.mean(matrix_ij)
    
        for i in range(M):
            for j in range(i, M):
                groupped_distance_matrix[j, i] = groupped_distance_matrix[i, j]
    
        return unique_labels, groupped_distance_matrix

    def plot_entropy_and_heatmap(self, save_prefix=None):
        
        unique_labels, adm = self.average_distance_matrix()
        linkage, order = hierarchical_clustering_order(adm)
        row_colors = list(map(lambda x: self.colormaps['core_annotation2']['colormap'].get(x), unique_labels[order]))

        g = sns.clustermap(-adm,
           row_colors=row_colors, cmap='RdYlBu_r')

        g.ax_heatmap.set_yticks(np.arange(adm.shape[0]) + 0.5)
        g.ax_heatmap.set_yticklabels([f'{lab} [{(self.adata.obs["core_annotation2"] == lab).sum()}]' for lab in unique_labels[order]], rotation=0)
        
        left, bottom, width, height = g.ax_heatmap.get_position().bounds
        
        # Add new axes on the right
        ax2 = g.fig.add_axes([left + width * 1.5, bottom, width * 0.7, height])  # [left, bottom, width, height]
        
        hp = sns.boxplot(data=self.adata.obs, x='entropy', y='core_annotation2', order=unique_labels[order],
             palette=self.colormaps['core_annotation2']['colormap'], ax=ax2)
        ax2.set_yticklabels(['' for _ in ax2.get_yticks()])
        ax2.set_ylabel('')
        # ax2.set_xticklabels(ax2.get_xticklabels(), rotation=90)
        if save_prefix is None:
            plt.show()
        else:
            plt.savefig(f"{save_prefix}.heatmap.pdf")

    def plot_tsne(self, columns=None, save_prefix=None):
        if columns is None:
            columns = ['extended_annotation2', 'major_component', 'leiden']
        
        for column in columns:
            fig, ax = plt.subplots(figsize=(5, 5))
            sc.pl.tsne(self.adata, color=column, show=False,
                       na_color=self.na_color, palette=self.colormaps.get(column, {'colormap': None})['colormap'], ax=ax)
            if save_prefix is None:
                plt.show()
            else:
                plt.savefig(f"{save_prefix}.{column}.pdf")
    
    def plot_mp(self, save_prefix=None):
        self.set_distances()
        self.calculate_clustering()
        self.set_colormaps()
        self.calculate_tsne()
        self.plot_entropy_and_heatmap(save_prefix)
        self.plot_tsne(save_prefix)

def hierarchical_clustering_order(matrix):
    Y = linkage(pdist(matrix, metric='euclidean'), method='single')
    return Y, leaves_list(Y)


def main(data, outpath):
    euclid_dist = pairwise_distances(data.obsm['embedding'].T, metric='jensenshannon')

    entropy = calc_entropy(euclid_dist, data.obs, 
        entropy_same_num=get_entropy_same_num(data.obs)
    )

    data.obs['entropy'] = entropy
    data.obsp['distance_matrix'] = euclid_dist

    mp = ModelPlotter(data, nearest_neighbors=50, resolution=2.5,
                    comp_prop=0.05)
    mp.plot_mp(save_prefix=outpath)


if __name__ == '__main__':
    samples_meta = pd.read_table(sys.argv[1])
    peaks_meta = pd.read_table(sys.argv[2])
    adata = ad.read_h5ad(sys.argv[3])
    mask123 = np.loadtxt(sys.argv[4], dtype=bool)

    adata.obs = samples_meta
    adata.obs['core_system'] = adata.obs['system']
    adata.obs['num_peaks_downsampled'] = adata.layers['binary'].sum(axis=1)

    adata.var = peaks_meta[mask123]
    main(adata, sys.argv[5])