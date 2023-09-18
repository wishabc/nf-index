import numpy as np
import pandas as pd
import sys
import json
from numba import jit
import perform_NMF as NMF
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
import subset_peaks
from sklearn.decomposition import PCA

class FeatureSelection:
    def __init__(self, params, signal_matrix, binary_matrix, sample_labels,
                 peaks_meta):
        self.sample_labels = sample_labels
        self.initial_signal_matrix = signal_matrix
        self.initial_binary_matrix = binary_matrix

        self.peaks_meta = peaks_meta

        self.params = params

        self.signal = self.binary = None
        self.mean_signal = self.mean_binary = None

    def select_peaks_for_clustering(self):
        self.calc_mean_matrices()
        self.set_confounders_mask()
        self.set_matrix_data()
        return self.add_peaks()

    def calc_mean_matrices(self):
        self.mean_signal, self.mean_binary = subset_peaks.average_matrices(
            self.initial_signal_matrix,
            self.initial_binary_matrix,
            self.sample_labels,
            reprod=self.params['Reproduced']
        )

    def set_confounders_mask(self):
        self.confounders_mask = self.peaks_meta.eval(
            self.params['Filtering_by_confounders']
        ) & (self.mean_binary.sum(axis=1) >= 1)
        mean_signal_filtered = self.initial_signal_matrix[self.confounders_mask, :].mean(axis=1)
        self.confounders_mask[self.confounders_mask] = mean_signal_filtered >= self.params["Filtering_by_mean_signal"]

    def set_matrix_data(self):
        if self.params['Calculate_gini_by'] == 'sample':
            self.signal = self.initial_signal_matrix[self.confounders_mask, :]
        elif self.params['Calculate_gini_by'] == 'group':
            self.signal = self.mean_signal[self.confounders_mask, :]
        else:
            raise ValueError

        if self.params['Add_peaks_by'] == 'sample':
            self.binary = self.initial_binary_matrix[self.confounders_mask, :]
        elif self.params['Add_peaks_by'] == 'group':
            self.binary = self.mean_binary[self.confounders_mask, :]
        else:
            raise ValueError

    def get_peaks_order(self):
        if self.params['Variance_metric'] == 'gini_index':
            means, gini, smoothed_gini_final, \
                gini_argsort, top_gini_mask = subset_peaks.get_gini_index_for_peaks(
                    self.signal, 1000)
            return gini_argsort
        elif self.params['Variance_metric'] == 'var':
            var = np.var(self.signal, axis=1)
            return np.argsort(var)[::-1]
        elif self.params['Variance_metric'] == 'explained_var':
            var = np.var(self.signal, axis=1)
            return np.argsort(var * self.peaks_meta['extended_annotation'][self.confounders_mask])[::-1]
        else:
            raise ValueError

    def add_peaks(self):
        sorted_peaks = self.get_peaks_order()
        new_mask_sub = subset_peaks.add_peaks(self.binary, sorted_peaks, 
                                            self.params['Add_peaks_mv'],
                                            self.params['Add_peaks_per_group'])
        mask = np.zeros(self.initial_signal_matrix.shape[0], dtype=bool)

        mask[self.confounders_mask] = new_mask_sub
        return mask


def minmax_norm(subset):
    return (subset - subset.min()) / (subset.max() - subset.min())


def pairwise_euclidean(X):
    return squareform(pdist(X.T, 'euclidean'))


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
            W, H, model = NMF.perform_NMF(self.data, n_components=self.method_params['n_components'])
            self.model_object = model
            self.embedding = W.T
            self.basis = H.T

    def plot_model(self):
        if self.method == 'pca':
            plt.plot(self.model_object.explained_variance_ratio_.cumsum()[:self.method_params['n_components']])
            plt.show()
            
            plt.plot(self.model_object.explained_variance_ratio_[:self.method_params['n_components']])
            plt.show()


@jit(nopython=True)
def accumulate_counts(sorted_labels, p_values, n):
    label_count = {}
    for i in range(n):
        label = sorted_labels[i]
        if label not in label_count:
            label_count[label] = 0.0
        label_count[label] += p_values[i]
    return label_count


def get_entropy_same_num(samples_meta):
    same = [{'Epithelial (Eye)', 'Epithelial'},
        {'Epithelial', 'Kidney'},
        {'Fetal brain', 'Neurovascular', 'Neuroblastoma/Retinoblastoma'},
        {'Heart', 'Fetal muscle'},
        {'Lung', 'Lung cancer cell line'},
        {'Fetal brain', 'Brain (Cerebellar)', 'Brain'},
        {'Pluripotent/Pluripotent-derived', 'Neuroblastoma/Retinoblastoma', 'Lung cancer cell line'},
        {'Myeloid progenitor', 'Myeloid leukemia', 'Myeloid leukemia (K562)'},
       {'T47D', 'MCF-7'},
       {'Kidney', 'Fetal kidney'},
       {'Lymphoblastoid cell line', 'T-cell', 'B-cell', 'NK cell'},
       {'M1 macrophage', 'M2 macrophage'}]
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


def pairwise_distances(X, metric='euclidean'):
    return squareform(pdist(X.T, metric))


def main(params, samples_meta, peaks_meta, signal_matrix, binary_matrix):
    sample_labels = np.unique(
        np.where(
            samples_meta['core_annotation2'].notna(),
            samples_meta['core_annotation2'],
            samples_meta['core_ontology_term']
        ), return_inverse=True
    )[1]
    fs = FeatureSelection(
        params=params,
        signal_matrix=signal_matrix,
        binary_matrix=binary_matrix,
        sample_labels=sample_labels,
        peaks_meta=peaks_meta
    )
    mask = fs.select_peaks_for_clustering()
    data, minv, maxv = minmax_norm(signal_matrix[mask, :])

    emb = Embedding(data, method='nmf', method_params={'n_components': params['N_components']})
    emb.set_default_params()
    emb.calculate_embedding()

    embedded_data = emb.embedding # / embedding.embedding.max(axis=1)[:, None]
    embedded_data = embedded_data / embedded_data.sum(axis=0)

    euclid_dist = pairwise_distances(embedded_data, metric='jensenshannon')

    entropy = calc_entropy(euclid_dist, samples_meta, entropy_same_num=get_entropy_same_num(samples_meta))

    return euclid_dist, entropy, mask, emb.embedding, emb.basis


if __name__ == '__main__':
    with open(sys.argv[1], 'r') as f:
        params = json.load(f)
    samples_meta = pd.read_table(sys.argv[2])
    peaks_meta = pd.read_table(sys.argv[3])

    assert 'frac_method' in peaks_meta.columns
    assert 'core_annotation2' in samples_meta.columns

    signal_matrix = np.load(sys.argv[4])
    binary_matrix = np.load(sys.argv[5])

    assert len(peaks_meta) == signal_matrix.shape[0] == binary_matrix.shape[0]
    assert len(samples_meta) == signal_matrix.shape[1] == binary_matrix.shape[1]

    euclid_dist, entropy, mask, W, H = main(params, samples_meta, peaks_meta, signal_matrix, binary_matrix)

    outdir = sys.argv[6]

    np.save(f'{outdir}.distances.npy', euclid_dist)
    np.save(f'{outdir}.entropy.npy', entropy)
    np.save(f'{outdir}.mask.npy', mask)
    np.save(f'{outdir}.H.npy', H)
    np.save(f'{outdir}.W.npy', W)
