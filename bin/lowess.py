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
                 seed_number=1832245,
                 sample_number=75_000,
                 bin_number=100,
                 min_peak_replication=0.25,
                 sample_method='raw',
                 cv_number=30,
                 jobs=1,
                 ):
        self.cv_number = cv_number
        self.min_peak_replication = min_peak_replication
        self.seed_number = seed_number
        self.sample_number = sample_number
        self.bin_number = bin_number
        self.peak_outlier_threshold = peak_outlier_threshold
        self.delta_fraction = delta_fraction
        self.correlation_limit = correlation_limit
        self.cv_fraction = self.delta = None
        self.sample_method = sample_method
        self.jobs = mp.cpu_count() if jobs == 0 else jobs
        self.set_randomizer()
    

    def set_randomizer(self):
        self.seed = np.random.RandomState(self.seed_number)
    
    def outlier_limit(self, x):
        """
        """
        if x.count() > self.sample_number:
            subset = x[self.sample_masked_array(x, size=self.sample_number)]
        else:
            subset = x
        subset = subset[~subset.mask]
        fitted = expon.fit(subset)

        return expon(fitted[0], fitted[1]).isf((1 - self.peak_outlier_threshold) / x.count())

    def sample_masked_array(self, arr, size):
        p = ~arr.mask
        return self.seed.choice(np.arange(arr.size)[p], size=int(size), replace=False)

    def select_peaks_uniform(self, peaks, decent_peaks_mask, ignore=None):
        """
        Returns row indices of selected peaks
        """
        if ignore is not None:
            peaks_mask = decent_peaks_mask & ~ignore
        else:
            peaks_mask = decent_peaks_mask

        masked_peaks = ma.masked_array(peaks, ~peaks_mask)
        k = min(self.sample_number, masked_peaks.count())

        if self.sample_method == 'random':
            result_indices = self.sample_masked_array(masked_peaks, k)
        else:
            if self.sample_method == 'log':
                vls = np.log(masked_peaks + 1)
            elif self.sample_method == 'raw':
                max_value = self.outlier_limit(masked_peaks)
                new_mask = ~masked_peaks.mask & (masked_peaks < max_value)
                vls = ma.masked_where(~new_mask, masked_peaks)
            else:
                raise ValueError('Method not in (random, log, raw)')
            bin_width = (vls.max() - vls.min()) / self.bin_number
            bin_size = np.ceil(k / self.bin_number)
            sampled_peaks_indicies = []

            for i in np.arange(vls.min(), vls.max(), bin_width):
                window_min = i
                window_max = i + bin_width
                new_mask = ~vls.mask & (vls >= window_min) & (vls < window_max)
                window_peaks = ma.masked_where(~new_mask, vls)

                if window_peaks.count() == 0:
                    continue
                sampled_window_peaks_indices = self.sample_masked_array(window_peaks,
                                                                        size=min(bin_size, window_peaks.count()))

                sampled_peaks_indicies.append(sampled_window_peaks_indices)

            result_indices = np.unique(np.concatenate(sampled_peaks_indicies))
        # Convert to mask
        res = np.zeros(peaks.shape, dtype=bool)
        res[result_indices] = True
        return res

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
            delta = self.delta_fraction * np.percentile(x, 99)

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


    def get_peak_subset(self, ref_peaks, num_samples_per_peak: np.ndarray,
                        density_mat: np.ndarray, weights: np.ndarray) -> np.ndarray:
        """
        Select a subset of peaks well correlated to a reference (mean or geometric mean)

        Returns:
        --------
        Indices for selected subset of peaks
        """
        n = np.max(num_samples_per_peak)
        perc = np.linspace(0, 1, 21)[:-1]
        i = np.where(perc >= self.min_peak_replication)[0]
        thresholds = np.floor(perc[i] * density_mat.shape[1])
        for ind, i in enumerate(thresholds):
            over = num_samples_per_peak >= i
            correlations = np.apply_along_axis(lambda x: spearmanr(x, ref_peaks[over])[0],
                axis=0,
                arr=density_mat[over, :]
            )
            avg_cor = np.average(correlations, weights=weights)
            logger.info(f'Selecting threshold, iteration #{ind}. Correlation {avg_cor}')
            if avg_cor > self.correlation_limit:
                break

        if i > n:
            logger.warning('Caution: individual samples may be poorly captured by mean!')

        return num_samples_per_peak >= i

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
    
    def get_xcounts(self, density_mat, mean_pseudocount, weights):
        logger.info('Computing mean and pseudocounts for each peak')
        mean_density = np.average(density_mat, axis=1, weights=weights)
        xvalues = np.log(mean_density + mean_pseudocount)
        return mean_density, xvalues

    def sample_peaks(self, density_mat: np.ndarray, mean_density: np.ndarray, peaks_mat: np.ndarray, weights: np.ndarray):
        """
        Select well-correlated peaks and sample a subset
        """
        num_samples_per_peak = self.get_num_samples_per_peak(peaks_mat)
        logger.info(f'Sampling representative (well-correlated) peaks (r2>{self.correlation_limit}) to mean')
        decent_peaks_mask = self.get_peak_subset(mean_density, num_samples_per_peak, density_mat, weights)
        sampled_peaks_mask = self.select_peaks_uniform(mean_density, decent_peaks_mask)
        logger.info(
            f'Found {decent_peaks_mask.sum():,} well-correlated peaks, using method "{self.sample_method}"'
            f' and sampled {sampled_peaks_mask.sum():,} peaks')

        self.delta = np.percentile(mean_density, 99) * self.delta_fraction

        return sampled_peaks_mask

    def fit_lowess_params(self, diffs: np.ndarray, xvalues: np.ndarray,
        sampled_peaks_mask: np.ndarray, weights: np.ndarray):
        _, S = diffs.shape
        logger.info('Computing LOWESS smoothing parameter via cross-validation')
        cv_set = self.seed.choice(S, size=min(self.cv_number, S), replace=False, p=weights)

        cv_fraction = np.mean(self.parallel_apply_2D(self.choose_fraction_cv, axis=0,
                                                     arr=diffs[:, cv_set], x=xvalues,
                                                     sampled=sampled_peaks_mask,
                                                     delta=self.delta))
        self.cv_fraction = cv_fraction

    def lowess_normalize(self, diffs: np.ndarray, xvalues: np.ndarray,
        sampled_peaks_mask: np.ndarray):
        """
        Normalizes to the mean of of the dataset
        Uses only well-correlated peaks to perform normalization
        """

        logger.info(f'Computing LOWESS on all the data with params - delta = {self.delta}, frac = {self.cv_fraction}')

        norm = self.parallel_apply_2D(self.fit_and_extrapolate, axis=0,
                                      arr=diffs, x=xvalues, sampled=sampled_peaks_mask,
                                      delta=self.delta, frac=self.cv_fraction)

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

    @staticmethod
    def get_scale_factors(matrix):
        return np.divide(1., (matrix.sum(axis=0) * 2 / 1.e5), dtype=dtype)


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
    
    dens_outpath = f'{base_path}.density.npy'
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

    weights = weights / weights.sum()
    data_norm = DataNormalize(jobs=p_args.jobs)
    scale_factors = data_norm.get_scale_factors(counts_matrix)
    density_matrix = counts_matrix * scale_factors
    pseudocounts = data_norm.get_pseudocounts(density_matrix)
    mean_pseudocount = pseudocounts.mean()
    np.save(dens_outpath, density_matrix)
    mean_scale_factors = scale_factors.mean()
    mat_and_pseudo = np.log(density_matrix + pseudocounts)
    del scale_factors
    del counts_matrix

    if model_params is None:
        mean_density, xvals = data_norm.get_xcounts(
            density_mat=density_matrix,
            mean_pseudocount=mean_pseudocount,
            weights=weights
        )
    
        sampled_mask = data_norm.sample_peaks(
            density_mat=density_matrix,
            mean_density=mean_density,
            peaks_mat=peaks_matrix,
            weights=weights)
        
        differences = (mat_and_pseudo.T - xvals).T
        del density_matrix
        del peaks_matrix
        gc.collect()

        data_norm.fit_lowess_params(diffs=differences,
            xvalues=xvals,
            sampled_peaks_mask=sampled_mask,
            weights=weights)
        deseq2_mean_sf = None
    else:
        del density_matrix
        del peaks_matrix
        gc.collect()

        xvals, sampled_mask, deseq2_mean_sf, weights = data_norm.load_params(model_params)
        differences = (mat_and_pseudo.T - xvals).T

    lowess_norm = data_norm.lowess_normalize(diffs=differences,
        xvalues=xvals, sampled_peaks_mask=sampled_mask)

    del differences
    gc.collect()
    np.save(lowess_outpath, lowess_norm)

    logger.info('Reconstructing normed matrix')
    density_matrix = np.load(dens_outpath)
    normed = density_matrix / lowess_norm
    normed = normed / mean_scale_factors

    del density_matrix
    del lowess_norm
    gc.collect()
    
    logger.info('Saving normed matrix')
    np.save(f'{base_path}.normed.npy', normed)
    logger.info('Reading raw tags...')
    counts_matrix = np.load(p_args.signal_matrix)

    deseq2_mean_sf = get_deseq2_scale_factors(counts_matrix,
        normed,
        f'{base_path}.scale_factors.npy',
        calculated_mean=deseq2_mean_sf,
        weights=weights
    )
    data_norm.save_params(model_save_params_path, xvals, sampled_mask, deseq2_mean_sf, weights)
