import numpy as np
import subset_peaks

class FeatureSelection:
    def __init__(self, params, signal_matrix, binary_matrix, sample_labels,
                 peaks_meta, signal_tr=1.5):
        self.sample_labels = sample_labels
        self.initial_signal_matrix = signal_matrix
        self.initial_binary_matrix = binary_matrix
        # Add asserts about analysis_df, signal_matrix, binary_matrix
        self.peaks_meta = peaks_meta
        self.signal_tr = signal_tr

        self.params = params

        self.signal = self.binary = None
        self.filtered_mean_signal = self.filtered_mean_binary = None

    def select_peaks_for_clustering(self):
        self.calc_mean_matrices()
        self.set_confounders_mask()
        self.set_matrix_data()
        return self.add_peaks()

    def calc_mean_matrices(self):
        self.filtered_mean_signal, self.filtered_mean_binary = subset_peaks.average_matrices(
            self.initial_signal_matrix,
            self.initial_binary_matrix,
            self.sample_labels,
            reprod=self.params['Reproduced']
        )

    def set_confounders_mask(self):
        self.confounders_mask = self.peaks_meta.eval(
            self.params['Filtering_by_confounders']
        ) & (self.filtered_mean_binary.sum(axis=1) >= 1)
        mean_signal_filtered = self.initial_signal_matrix[self.confounders_mask, :].mean(axis=1)
        self.confounders_mask[self.confounders_mask] = mean_signal_filtered >= self.signal_tr


    def set_matrix_data(self):
        if self.params['Calculate_gini_by'] == 'sample':
            self.signal = self.initial_signal_matrix[self.confounders_mask, :]
        elif self.params['Calculate_gini_by'] == 'group':
            self.signal = self.filtered_mean_signal
        else:
            raise ValueError

        if self.params['Add_peaks_by'] == 'sample':
            self.binary = self.initial_binary_matrix[self.confounders_mask, :]
        elif self.params['Add_peaks_by'] == 'group':
            self.binary = self.filtered_mean_binary
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
        new_mask_sub = subset_peaks.add_peaks(self.binary, sorted_peaks, self.params['Add_peaks_mv'],
                                              self.params['Add_peaks_per_group'])
        mask = np.zeros(self.initial_signal_matrix.shape[0], dtype=bool)
        mask[self.confounders_mask] = new_mask_sub
        return mask


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
    data = minmax_norm(signal_matrix[mask, :])
    euclid_dist = pairwise_euclidean(data)
    entropy = calc_entropy(euclid_dist, samples_meta, entropy_same_num=)
    euclid_dist ## saveme
    entropy ## saveme
    mask ## saveme

def calc_entropy(euclid_dist, samples_meta, entropy_same_num):

    return get_entropy_scores(
        euclid_dist, 
        np.unique(samples_meta['core_annotation2'], return_inverse=True)[1],
        entropy_same_num
    )


def plot_model(analysis_df, colors_order, label_counts, n_peaks):
    mean = np.mean(analysis_df.groupby('core_annotation2')['delta_entropy'].apply(np.mean))
    median = np.median(analysis_df.groupby('core_annotation2')['delta_entropy'].apply(np.median))
    fig, ax = plt.subplots(figsize=(12, 4))
    hp = sns.boxplot(data=analysis_df, y='delta_entropy', x='core_annotation2',
                     order=label_counts.index, palette=colors_order)
    ax.axhline(0, ls='--', color='grey', lw=2)
    ax.set_xticklabels(
        [f'{row["core_annotation2"]} [{row["count"]}]' for _, row in label_counts.reset_index().iterrows()],
        rotation=90)
    plt.title(f'{n_peaks} peaks, net_effect {median:.3f}, {mean:.3f}')
    plt.show()