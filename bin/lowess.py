import os
import numpy as np
from numpy import ma
from scipy.stats import expon, spearmanr
from statsmodels.nonparametric import smoothers_lowess
from convert_to_numpy import read_matrix
import multiprocessing as mp
import logging
import sys
import argparse
import gc
import json

handler = logging.StreamHandler(sys.stdout)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s  %(levelname)s  %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

dtype = np.float64


class DataNormalize:
    def __init__(self,
                 peak_outlier_threshold=0.999,
                 delta_fraction=0.001,
                 correlation_limit=0.8,
                 seed=42,
                 bin_number=100,
                 min_peak_replication=0.25,
                 sample_method='log',
                 cv_number=30,
                 common_scale=50_000,
                 jobs=1,
                 ):
        self.cv_number = cv_number
        self.min_peak_replication = min_peak_replication
        self.seed_number = seed
        self.bin_number = bin_number
        self.peak_outlier_threshold = peak_outlier_threshold
        self.delta_fraction = delta_fraction
        self.correlation_limit = correlation_limit
        self.cv_fraction = None
        self.delta = 0
        self.sample_method = sample_method
        self.jobs = mp.cpu_count() if jobs == 0 else jobs
        self.common_scale = common_scale
        self.set_randomizer()
    

    def set_randomizer(self):
        self.seed = np.random.RandomState(self.seed_number)
    
    def outlier_limit(self, x):
        """
        Returns the outlier limit for the data by fitting an exponential distribution to the right tail
        of the data (99% quantile)
        """
        fit_q = 0.8
        x = x[~x.mask]
        if self.peak_outlier_threshold < fit_q:
            return np.quantile(x, self.peak_outlier_threshold)
        elif self.peak_outlier_threshold == 1:
            return np.inf
        fit_q_val = np.quantile(x, fit_q)
        values_to_fit = np.exp(x)[x >= fit_q_val]
        fitted = expon.fit(values_to_fit, floc=values_to_fit.min())
        return np.log(expon(*fitted).isf(1 - (self.peak_outlier_threshold - fit_q) / (1 - fit_q)))

    def sample_masked_array(self, arr, size):
        p = ~arr.mask
        return self.seed.choice(np.arange(arr.size)[p], size=int(size), replace=False)

    @staticmethod
    def check_argsort_thresholds(argsorts, bin_size):
        unique_elements, counts = np.unique(argsorts, return_counts=True)
        max_elems_in_bin = unique_elements[counts == counts.max()].max()

        if max_elems_in_bin < bin_size:
            logger.warning(f'Not enough peaks to sample in one of the bins: {max_elems_in_bin} < {bin_size}')
            #bin_size = (max_elems_in_bin + bin_size) / 2

    @staticmethod
    def masked_ranks(a):
        ranks = np.ma.array(np.empty_like(a, dtype=int), mask=a.mask)
        ranks[~a.mask] = np.argsort(np.argsort(a[~a.mask]))
        return ranks / a.count()

    @staticmethod
    def weighted_variance(x, w):
        """
        Returns weighted variance
        """
        mean = np.average(x, weights=w, axis=1)
        return np.average((x - mean) ** 2, weights=w, axis=1)

    def select_peaks_uniform(self, log_cpm, mean_log_cpm, weights, num_samples_per_peak):
        """
        Returns row indices of selected peaks
        """
        reproducible_peaks = num_samples_per_peak >= self.min_peak_replication * log_cpm.shape[1]
    
        repr_log_means = ma.masked_array(mean_log_cpm, ~reproducible_peaks)

        max_value = self.outlier_limit(repr_log_means)
        masked_log_means = ma.masked_where(
            ~repr_log_means.mask & (repr_log_means < max_value), 
            repr_log_means
        )

        bin_edges = np.linspace(masked_log_means.min(), masked_log_means.max(), self.bin_number + 1)
        
        sampled_peaks_indicies = np.zeros(mean_log_cpm.shape, dtype=bool)
        
        peak_variance = self.weighted_variance(log_cpm)
        
        per_bin_ranks = np.full_like(masked_log_means, np.nan)
        per_bin_ranks = np.ma.masked_array(per_bin_ranks, masked_log_means.mask)
        for i in range(self.bin_number):
            window_min = bin_edges[i]
            window_max = bin_edges[i + 1]

            new_mask = ~masked_log_means.mask & (masked_log_means >= window_min) & (masked_log_means < window_max)
            if i == self.bin_number - 1:
                new_mask |= (masked_log_means == window_max)
            peak_variance_window = ma.masked_where(~new_mask, peak_variance)
            per_bin_ranks[new_mask] = self.masked_ranks(peak_variance_window)[new_mask]

        #self.check_argsort_thresholds(per_bin_argsorts, bin_size)
    
        # top_by_variance_thresholds = np.linspace(0, bin_size, 10)[::-1]
        bottom_by_variance_thresholds = np.linspace(0, 1, 20)[1:][::-1]
        sampled_peaks_indicies = self.choose_best_score_by_correlation(
            mean_log_cpm=masked_log_means,
            log_cpm=log_cpm,
            peak_scores=-per_bin_ranks, # - to choose "top" peaks
            score_thresholds=-bottom_by_variance_thresholds,
            weights=weights,
        )
        return sampled_peaks_indicies

    @staticmethod
    def extrapolate(interpolated, to_extrap, base, predict):
        """
        Sets the values out of the range to the edge values.
        For logged data, equivalent to linear approximation between 0 and the point
        """

        sample_order = np.argsort(base)
        min_value = np.min(base)
        max_value = np.max(base)

        under_min = to_extrap < min_value
        over_max = to_extrap > max_value

        min_predict = predict[sample_order[0]]
        max_predict = predict[sample_order[-1]]

        interpolated[under_min] = min_predict
        interpolated[over_max] = max_predict

        return interpolated

    def get_extrapolation(self, x, lowess_est, sampled):
        """
        """
        sort_ind = np.argsort(x[sampled])
        interpolated = np.interp(x, x[sampled][sort_ind], lowess_est[sort_ind])
        extrapolated = self.extrapolate(interpolated, x, x[sampled], lowess_est)
        return extrapolated

    def run_lowess(self, y, x, sampled, frac, delta):
        return smoothers_lowess.lowess(y[sampled], x[sampled],
                                                return_sorted=False, it=4,
                                                frac=frac, delta=delta)

    def fit_and_extrapolate(self, y, x, sampled, frac, delta):
        smoothed_values = self.run_lowess(y, x, sampled, frac, delta)
        return self.get_extrapolation(x, smoothed_values, sampled)

    def choose_fraction_cv(self, y, x, sampled, start=0.1, end=0.8, step=0.1, delta=None):
        """
        choose smoothing parameter for lowess by cross-validation
        sampled:
        """

        lo = x[sampled].min()
        hi = x[sampled].max()

        within_range_indices = (x >= lo) & (x <= hi)
        cv_sample = within_range_indices & ~sampled

        min_err = np.inf
        best_frac = 0

        if delta is None:
            delta = self.delta_fraction * (np.percentile(x, 99) - np.percentile(x, 1))

        for frac in np.arange(start, end + step, step):
            interpolated = self.fit_and_extrapolate(y, x, sampled, frac, delta)

            if np.isnan(interpolated.max()):
                err = np.inf
            else:
                err = np.square(interpolated[cv_sample] - y[cv_sample]).mean()

            if err < min_err:
                min_err = err
                best_frac = frac

        return best_frac

    @staticmethod
    def get_num_samples_per_peak(matrix):
        """
        Returns number of samples for each peak
        """
        return matrix.sum(axis=1)

    @staticmethod
    def get_pseudocounts(matrix):
        """
        Compute pseudocounts for each sample in the matrix
        """
        return np.nanmin(ma.masked_where(matrix <= 0.0, matrix), axis=0)

    def choose_best_score_by_correlation(self, 
                        mean_log_cpm: np.ma.MaskedArray,
                        log_cpm: np.ndarray,
                        peak_scores: np.ma.MaskedArray,
                        score_thresholds: np.ndarray,
                        weights: np.ndarray) -> np.ndarray:
        """
        Select a subset of peaks well correlated to a reference (mean or geometric mean)

        Returns:
        --------
        Indices for selected subset of peaks
        """
        best_thr = score_thresholds[0]
        best_correlation = -np.inf
        for ind, thr in enumerate(score_thresholds):
            over = (peak_scores >= thr).filled(False)
            correlations = np.apply_along_axis(
                lambda x: spearmanr(x, mean_log_cpm[over])[0],
                axis=0,
                arr=log_cpm[over, :]
            )
            avg_cor = np.average(correlations, weights=weights)
            if avg_cor > best_correlation:
                best_correlation = avg_cor
                best_thr = thr
            logger.info(f'Selecting threshold, iteration #{ind}. Correlation {avg_cor}.')
            if avg_cor > self.correlation_limit:
                break
        else:
            logger.warning(f'Caution: individual samples may be poorly captured by mean! Best correlation({best_correlation:.2f})')

        return peak_scores >= best_thr

    @staticmethod
    def unpack_and_apply(x, func1d, kwargs, *args):
        return func1d(x, *args, **kwargs)

    def parallel_apply_2D(self, func1d, axis, arr, *args, **kwargs):
        """
        Parallel version of apply_along_axis() for 2D matrices
        """
        if axis == 1:
            arr = arr.T
        n_samples = arr.shape[1]
        jobs = min(self.jobs, n_samples)
        if jobs > 1:
            ctx = mp.get_context("forkserver")
            with ctx.Pool(jobs) as p:
                individual_results = p.starmap(self.unpack_and_apply,
                                               [(arr[:, index], func1d, kwargs, *args) for
                                                index in range(n_samples)])

            if np.ndim(individual_results) > 1:
                result = np.vstack(individual_results).T
            else:
                result = np.array(individual_results)
            return result if axis == 0 else result.T
        else:
            return np.apply_along_axis(func1d, axis, arr, *args, **kwargs)

    def sample_peaks(self, log_cpm: np.ndarray, mean_log_cpm: np.ndarray, num_samples_per_peak: np.ndarray, weights: np.ndarray):
        """
        Select well-correlated peaks and sample a subset
        """
        logger.info(f'Sampling representative (well-correlated) peaks (r2>{self.correlation_limit}) to mean')
    
        sampled_peaks_mask = self.select_peaks_uniform(
            log_cpm, 
            mean_log_cpm, 
            weights=weights,
            num_samples_per_peak=num_samples_per_peak
        )
        logger.info(f'Sampled {sampled_peaks_mask.sum():,} well-correlated peaks')

        return sampled_peaks_mask

    def fit_lowess_params(self, diffs: np.ndarray, xvalues: np.ndarray, sampled_peaks_mask: np.ndarray, weights: np.ndarray):
        _, S = diffs.shape
        logger.info('Computing LOWESS smoothing parameter via cross-validation')
        cv_set = self.seed.choice(S, size=min(self.cv_number, S), replace=False, p=weights)
        
        cv_fraction = np.mean(
            self.parallel_apply_2D(
                self.choose_fraction_cv,
                axis=0,
                arr=diffs[:, cv_set], 
                x=xvalues,
                sampled=sampled_peaks_mask,
                delta=self.delta
            )
        )
        self.cv_fraction = cv_fraction

    def lowess_normalize(self, diffs: np.ndarray, xvalues: np.ndarray,
        sampled_peaks_mask: np.ndarray):
        """
        Normalizes to the mean of of the dataset
        Uses only well-correlated peaks to perform normalization
        """
        logger.info(f'Computing LOWESS on all the data with params - delta = {self.delta}, frac = {self.cv_fraction}')

        norm = self.parallel_apply_2D(
            self.fit_and_extrapolate, 
            axis=0,
            arr=diffs,
            x=xvalues,
            sampled=sampled_peaks_mask,
            delta=self.delta,
            frac=self.cv_fraction
        )

        logger.info('Normalizing finished')
        return np.exp(norm)

    def load_params(self, model_params):
        # unpack params from npz array
        with open(f'{model_params}.json') as f:
            class_params = json.load(f)
        for key, value in class_params.items():
            setattr(self, key, value)
        self.set_randomizer()
        arrays = np.load(f'{model_params}.npz')
        return arrays['xvalues'], arrays['sampled_mask'], arrays['deseq2_mean_sf'], arrays['weights']
        
    def save_params(self, save_path, xvals, sampled_mask, deseq2_mean_sf, weights):
        for ext in '.npz', '.json':
            if os.path.exists(f'{save_path}{ext}'):
                logger.warning(f'File {save_path}{ext} exists, model params were not saved')
                return
        
        with open(f'{save_path}.json', 'w') as f:
            params_dict  = {x:y for x,y in self.__dict__.items() if x != 'seed'}
            json.dump(params_dict, f, indent=2)

        np.savez_compressed(f'{save_path}.npz',
            xvalues=xvals,
            sampled_mask=sampled_mask,
            deseq2_mean_sf=deseq2_mean_sf,
            weights=weights
        )

    def get_scale_factors(self, matrix):
        # Divide by 50k - average # of peaks per sample
        return np.divide(self.common_scale, matrix.sum(axis=0), dtype=dtype)


def check_and_open_matrix_file(path, outpath):
    if not os.path.exists(path):
        raise ValueError(f'{path} do not exist')
    _, ext = os.path.splitext(path)
    if ext == '.npy':
        return np.load(path)
    else:
        np_arr = read_matrix(path)
        if np.ma.isMA(np_arr):
            np_arr = np_arr.filled(0)
        np.save(outpath, np_arr)
        return np_arr



def get_deseq2_scale_factors(raw_tags, as_normed_matrix, scale_factor_path, calculated_mean, weights):
    logger.info('Calculating scale factors...')
    sf = raw_tags / as_normed_matrix
    i = (raw_tags == 0) | (as_normed_matrix == 0)
    sf[i] = 1

    sf[np.isnan(sf)] = 1
    sf[~np.isfinite(sf)] = 1
    if calculated_mean is not None:
        sf_geomean = calculated_mean
    else:
        sf_geomean = np.exp(np.average(np.log(sf), axis=1, weights=weights))
    as_scale_factors = (sf.T / sf_geomean).T
    np.save(scale_factor_path, as_scale_factors)
    return sf_geomean


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Matrix normalization using lowess')
    parser.add_argument('peak_matrix', help='Path to binary peaks matrix')
    parser.add_argument('signal_matrix', help='Path to matrix with read counts for each peak in every sample')
    parser.add_argument('output', help='Path to directory to save normalized matrix into.')
    parser.add_argument('--prefix', help='Filenames prefix', default='matrix')
    parser.add_argument('--weights', help='Path to weights (for weighted lowess)')
    parser.add_argument('--jobs', type=int,
                        help='Number of jobs to parallelize calculations '
                             '(can\'t be larger than number of samples. If 0 is provided - uses all available cores')
    parser.add_argument('--model_params', help='Use existing lowess params ffor normalization. Expected to provide basename of the files with lowess params. E.g. provide "prefix.lowess_params" to load both params files: "prefix.lowess_params.json" and "prefix.lowess_params.npz".', default=None)
    p_args = parser.parse_args()

    model_params = p_args.model_params

    if not os.path.exists(p_args.output):
        os.mkdir(p_args.output)
    base_path = os.path.join(p_args.output, p_args.prefix)
    
    cpm_matrix_outpath = f'{base_path}.density.npy'
    lowess_outpath = f'{base_path}.lowess.npy'
    
    model_save_params_path = f'{base_path}.lowess_params'

    logger.info('Reading matrices')
    counts_matrix = np.load(p_args.signal_matrix)
    peaks_matrix = np.load(p_args.peak_matrix)
    
    N, S = counts_matrix.shape
    assert counts_matrix.shape == peaks_matrix.shape
    logger.info(f'Normalizing matrix with shape: {N:,};{S}')

    if p_args.weights is not None:
        weights = np.load(p_args.weights)
    else:
        weights = np.ones(S)

    pseudocount = 1
    weights = weights / weights.sum()
    data_norm = DataNormalize(jobs=p_args.jobs)
    scale_factors = data_norm.get_scale_factors(counts_matrix)
    cpm_matrix = counts_matrix * scale_factors
    np.save(cpm_matrix_outpath, cpm_matrix)
    log_cpm_matrix = np.log(cpm_matrix + pseudocount * scale_factors)

    num_samples_per_peak = data_norm.get_num_samples_per_peak(peaks_matrix)
    del scale_factors
    del counts_matrix
    del peaks_matrix

    if model_params is None:
        mean_log_cpm = np.average(log_cpm_matrix, axis=1, weights=weights)

        sampled_mask = data_norm.sample_peaks(
            log_cpm=log_cpm_matrix,
            mean_log_cpm=mean_log_cpm,
            num_samples_per_peak=num_samples_per_peak,
            weights=weights
        )
        
        log_differences = log_cpm_matrix - mean_log_cpm[:, None]
        
        data_norm.fit_lowess_params(
            diffs=log_differences,
            xvalues=mean_log_cpm,
            sampled_peaks_mask=sampled_mask,
            weights=weights
        )
        deseq2_mean_sf = None
    else:
        mean_log_cpm, sampled_mask, deseq2_mean_sf, weights = data_norm.load_params(model_params)
        log_differences = log_cpm_matrix - mean_log_cpm[:, None]

    del log_cpm_matrix
    gc.collect()
    lowess_norm = data_norm.lowess_normalize(
        diffs=log_differences,
        xvalues=mean_log_cpm,
        sampled_peaks_mask=sampled_mask
    )

    del log_differences
    gc.collect()
    np.save(lowess_outpath, lowess_norm)

    logger.info('Reconstructing normed matrix')
    cpm_matrix = np.load(cpm_matrix_outpath)
    normed = cpm_matrix / lowess_norm
    # Not sure
    normed = normed / data_norm.common_scale

    del cpm_matrix
    del lowess_norm
    gc.collect()
    
    logger.info('Saving normed matrix')
    np.save(f'{base_path}.normed.npy', normed)
    logger.info('Reading raw tags...')
    counts_matrix = np.load(p_args.signal_matrix)

    deseq2_mean_sf = get_deseq2_scale_factors(
        counts_matrix,
        normed,
        f'{base_path}.scale_factors.npy',
        calculated_mean=deseq2_mean_sf,
        weights=weights
    )
    data_norm.save_params(
        model_save_params_path, 
        mean_log_cpm,
        sampled_mask,
        deseq2_mean_sf,
        weights
    )
