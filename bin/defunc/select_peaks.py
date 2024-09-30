import numpy as np
import pandas as pd
import sys
import json
from numba import jit
import perform_NMF as pNMF
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
import anndata as ad


def average_matrices(normalized_matrix, binary_matrix, meta_labels, reprod=1):
    unique_labels = np.unique(meta_labels)
    new_norm_matrix = np.zeros(shape=(normalized_matrix.shape[0], unique_labels.shape[0]),
                                dtype=normalized_matrix.dtype)
    new_binary_matrix = np.zeros(shape=(binary_matrix.shape[0], unique_labels.shape[0]),
                                dtype=bool)
    for label in unique_labels:
        new_norm_matrix[:, label] = normalized_matrix[:, meta_labels == label].mean(axis=1)
        new_binary_matrix[:, label] = binary_matrix[:, meta_labels == label].sum(axis=1) >= np.minimum(reprod, (meta_labels == label).sum())

    return new_norm_matrix, new_binary_matrix


def add_peaks(binary_matrix, gini_argsort, num_peaks=0, min_peaks_per_sample=0):
    sorted_binary = binary_matrix[gini_argsort, :]
    sorted_binary = np.where(sorted_binary.cumsum(axis=0) > min_peaks_per_sample, False, sorted_binary)
    sorted_binary[:num_peaks, :] = True

    inv_argsort = np.zeros_like(gini_argsort)
    inv_argsort[gini_argsort] = np.arange(len(gini_argsort))

    return sorted_binary.any(axis=1)[inv_argsort]


def add_sample_labels(adata, by='sample_label', fallback='core_ontology_term'):
    adata.obs['group'] = np.unique(
        np.where(
            adata.obs[by].notna(),
            adata.obs[by],
            'NA' if fallback == 'NA' else adata.obs[fallback]
        ), return_inverse=True
    )[1]

def calc_mean_matrices(adata, rep):
    adata.varm['averaged_signal_matrix'], adata.varm['averaged_binary_matrix'] = average_matrices(
        adata.X.T,
        adata.layers['binary'].T,
        adata.obs['group'],
        reprod=rep
    )
    adata.uns['Reproduced'] = rep


class Embedding():
    def __init__(self, data, method='raw', method_params={}, random_state=234234):
        if method not in ['raw', 'pca', 'nmf']:
            raise NotImplementedError
        self.method = method
        self.data = data
        self.method_params = method_params
        self.random_state = random_state

        self.embedding = None
        self.model_object = None

    def set_default_params(self):
        if self.method == 'pca':
            self.method_params.setdefault('n_components', 80)

    def calculate_embedding(self):
        if self.method == 'raw':
            self.embedding = self.data
        elif self.method == 'pca':
            self.model_object = PCA(
                n_components=self.method_params['n_components'],
                random_state=self.random_state
            )
            self.embedding = self.model_object.fit_transform(self.data.T, ).T
        elif self.method == 'nmf':
            W, H, model = pNMF.perform_NMF(self.data, n_components=self.method_params['n_components'])
            self.model_object = model
            self.embedding = W.T

    def plot_model(self):
        if self.method == 'pca':
            plt.plot(self.model_object.explained_variance_ratio_.cumsum()[:self.method_params['n_components']])
            plt.show()
            
            plt.plot(self.model_object.explained_variance_ratio_[:self.method_params['n_components']])
            plt.show()


class FeatureSelection:
    def __init__(self, params, adata):
        self.adata = adata
        self.params = params

        self.signal = self.binary = None

        self.adata.var['variance'] = None

        self.variance = None

    def select_peaks_for_clustering(self):
        # if (self.adata.uns['Reproduced'] != self.params['Reproduced'] or
        #     'averaged_signal_matrix' not in self.adata.varm or
        #     'averaged_binary_matrix' not in self.adata.varm):
        calc_mean_matrices(self.adata, self.params['Reproduced'])
        self.set_confounders_mask()
        self.set_matrix_data()
        self.calculate_variance()
        self.set_peak_ranks()
        return self.add_peaks()

    def set_confounders_mask(self):
        confounders_mask = self.adata.var.eval(
            self.params['Filtering_by_confounders']
        ) & (self.adata.varm['averaged_binary_matrix'].sum(axis=1) >= 1)
        mean_density = self.adata.X.T.mean(axis=1)
        mean_signal_filtered = mean_density[confounders_mask]
        thr = get_threshold(mean_signal_filtered, self.params["Filtering_by_mean_signal"])
        confounders_mask[confounders_mask] = mean_signal_filtered >= thr
        
        self.adata.var['confounders_mask'] = confounders_mask
        self.adata.var['mean_density'] = mean_density

        self.filtered_adata = self.adata[:, self.adata.var['confounders_mask']]

    def set_matrix_data(self):
        if self.params['Calculate_gini_by'] == 'sample':
            self.signal = self.filtered_adata.X.T
        elif self.params['Calculate_gini_by'] == 'group':
            self.signal = self.filtered_adata.varm['averaged_signal_matrix']
        else:
            raise ValueError

        if self.params['Add_peaks_by'] == 'sample':
            self.binary = self.filtered_adata.layers['binary'].T
        elif self.params['Add_peaks_by'] == 'group':
            self.binary = self.filtered_adata.varm['averaged_binary_matrix']
        else:
            raise ValueError

    def calculate_variance(self):
        print('Calculating variance')
        if self.params['Variance_metric'] in ['var', 'binned_var']:
            variance = np.var(self.signal, axis=1)
        elif self.params['Variance_metric'] == 'explained_var':
            variance = np.var(self.signal, axis=1) * self.filtered_adata.var['sample_label']
        else:
            raise ValueError

        with pd.option_context('mode.chained_assignment', None):
            self.variance = self.adata.var['variance'][self.adata.var['confounders_mask']] = variance

    def set_peak_ranks(self):
        if self.params['Variance_metric'] == 'binned_var':
            self.peak_ranks = self.get_peak_ranks_binned()
        else:
            self.peak_ranks = np.argsort(self.variance)[::-1]

    def get_peak_ranks_binned(self):
        ranks = np.zeros(self.variance.shape, dtype=bool)

        bins = np.quantile(self.filtered_adata.var['mean_density'], np.linspace(0, 1, 11))
        df = pd.DataFrame({'mean_signal_bin': pd.cut(self.filtered_adata.var['mean_density'], bins),
                     'variance': self.variance})
        for bin_label in df['mean_signal_bin'].cat.categories:
            bin_idx = df['mean_signal_bin'] == bin_label
            ranks[bin_idx] = np.argsort(df[bin_idx]['variance']) / len(bin_idx)
        
        return np.argsort(ranks)[::-1]

    def add_peaks(self):
        print('Adding peaks')
        num_peaks_to_leave = int(self.params['Add_peaks_mv'] * len(self.peak_ranks))
        new_mask_sub = add_peaks(self.binary, self.peak_ranks, num_peaks_to_leave,
                                              self.params['Add_peaks_per_group'])
        mask = np.zeros(self.adata.shape[1], dtype=bool)
        mask[self.adata.var['confounders_mask']] = new_mask_sub
        self.adata.var['most_variable_peaks'] = mask
        return mask

def minmax_norm(subset, minv=None, maxv=None):
    if minv is None:
        minv = subset.min()
    if maxv is None:
        maxv = subset.max()
    return (subset - minv) / (maxv - minv), minv, maxv


# def pairwise_distances(X, metric='euclidean'):
#     return squareform(pdist(X.T, metric))


# def calc_entropy(euclid_dist, samples_meta, entropy_same_num):
#     validation_idx = samples_meta['core_annotation'].notna()

#     entropy = np.full(validation_idx.shape, np.nan)
#     entropy[validation_idx] = get_entropy_scores(
#         euclid_dist[validation_idx, :][:, validation_idx], 
#         np.unique(samples_meta.loc[validation_idx, 'core_annotation'], return_inverse=True)[1],
#         entropy_same_num
#     )

#     return entropy


# def get_entropy_scores(euclid, labels_array, same_num_dict=None):
#     annotation_entropy_scores = []
#     for i in range(len(labels_array)):
#         annotation_entropy_scores.append(get_metric(euclid[i, :], labels_array, same_num_dict=same_num_dict))
    
#     annotation_entropy_scores = np.array(annotation_entropy_scores)
    
#     return annotation_entropy_scores


# def get_metric(distances, labels, same_num_dict=None):
#     n = len(distances)
#     assert len(labels) == n
    
#     sorted_indices = np.argsort(distances)
#     sorted_labels = labels[sorted_indices]

#     n_nearest = (sorted_labels == sorted_labels[0]).sum() * 0.5
#     p = np.power(0.1, 1/n_nearest)
    
#     p_values = np.power(p, np.arange(n))
#     p_values /= p_values.sum()

#     # Ignore similar annotations in entropy score, consider them same as original
#     if same_num_dict is not None and sorted_labels[0] in same_num_dict.keys():
#         for v in same_num_dict[sorted_labels[0]]:
#             sorted_labels[sorted_labels == v] = sorted_labels[0]
    
    
#     label_count = accumulate_counts(sorted_labels, p_values, n)
    
#     values = np.array(list(label_count.values()))
#     values = values[values != 0]
#     entropy = -np.sum(values * np.log2(values))
    
#     return entropy


def matrices_to_adata(signal_matrix, binary_matrix, samples_meta, peaks_meta):
    assert binary_matrix.shape[0] == peaks_meta.shape[0]
    assert binary_matrix.shape[1] == samples_meta.shape[0]

    assert signal_matrix.shape[0] == peaks_meta.shape[0]
    assert signal_matrix.shape[1] == samples_meta.shape[0]
    
    adata = ad.AnnData(signal_matrix.T)
    adata.obs_names = samples_meta.index
    adata.var_names = peaks_meta.index
    adata.layers['binary'] = binary_matrix.T
    adata.obs = samples_meta
    adata.var = peaks_meta
    return adata


@jit(nopython=True)
def accumulate_counts(sorted_labels, p_values, n):
    label_count = {}
    for i in range(n):
        label = sorted_labels[i]
        if label not in label_count:
            label_count[label] = 0.0
        label_count[label] += p_values[i]
    return label_count


def get_cumul_fraction(array, normalize=True):
    x = np.sort(array)[::-1]
    y = np.arange(len(x)) / (float(len(x)) if normalize else 1_000)
    return x, y

def get_threshold(array, thr):
    x, y = get_cumul_fraction(array)
    idx = min(np.searchsorted(y, thr), len(y) - 1)
    return x[idx]


def main(params, samples_meta, peaks_meta, signal_matrix, binary_matrix):
    adata = matrices_to_adata(signal_matrix, binary_matrix, samples_meta, peaks_meta)
    add_sample_labels(adata, by='core_annotation')
    calc_mean_matrices(adata, rep=1) 
    print('Mean matrices calculated')  
    fs = FeatureSelection(
        params=params,
        adata=adata
    )
    mask = fs.select_peaks_for_clustering()
    print('Peaks selected')
    filtered_adata = adata[:, mask]
    filtered_adata.layers['minmax'], minv, maxv = minmax_norm(filtered_adata.X)
    filtered_adata.uns['minv'] = minv
    filtered_adata.uns['maxv'] = maxv

    emb_handler = Embedding(filtered_adata.layers['minmax' if params['Matrix'] == 'signal' else 'binary'].astype(float).T,
        method='nmf',
        method_params={'n_components': params['N_components']})
    emb_handler.set_default_params()
    # emb_handler.calculate_embedding()

    # embedding = emb_handler.embedding

    # filtered_adata.obsm['embedding'] = embedding.T
    # filtered_adata.obsm['normalized_embedding'] = (embedding / embedding.sum(axis=0)).T

    # filtered_adata.obsp['distance_matrix'] = pairwise_distances(
    #     filtered_adata.obsm['normalized_embedding'].T,
    #     metric='jensenshannon'
    # )
    
    # filtered_adata.varm['basis'] = emb_handler.model_object.components_.T
    return filtered_adata, mask


if __name__ == '__main__':
    with open(sys.argv[1], 'r') as f:
        params = json.load(f)
    samples_meta = pd.read_table(sys.argv[2])
    peaks_meta = pd.read_table(sys.argv[3])

    signal_matrix = np.load(sys.argv[4])
    binary_matrix = np.load(sys.argv[5])

    assert len(peaks_meta) == signal_matrix.shape[0] == binary_matrix.shape[0]
    assert len(samples_meta) == signal_matrix.shape[1] == binary_matrix.shape[1]
    print('Matrices loaded')
    adata, mask = main(params, samples_meta, peaks_meta, signal_matrix, binary_matrix)

    adata.obs.drop(columns=adata.obs.columns, inplace=True)
    adata.var.drop(columns=adata.var.columns, inplace=True)
    adata.write(f"{sys.argv[6]}.h5")
    print('Data saved')
    np.savetxt(f"{sys.argv[6]}.selected_peaks.mask.txt", mask)
    
