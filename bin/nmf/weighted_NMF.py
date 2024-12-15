import time
import numpy as np
from sklearn.decomposition import NMF
from sklearn.decomposition._nmf import _beta_loss_to_float, trace_dot
import numpy as np
import scipy.sparse as sp
import time
import warnings


from sklearn._config import config_context
from sklearn.exceptions import ConvergenceWarning
from sklearn.utils.extmath import safe_sparse_dot, squared_norm
from sklearn.utils.validation import check_is_fitted, check_non_negative


EPSILON = np.finfo(np.float32).eps


class WeightedNMF(NMF):

    def fit_transform(self, X, y=None, W=None, H=None, W_weights=None, H_weights=None):

        self._validate_params()

        X = self._validate_data(
            X, accept_sparse=("csr", "csc"), dtype=[np.float64, np.float32]
        )
        if W_weights is not None and self.beta_loss != 'frobenius':
            raise NotImplementedError
        with config_context(assume_finite=True):
            W, H, n_iter = self._fit_transform(X, W=W, H=H, W_weights=W_weights, H_weights=H_weights)

        self.reconstruction_err_ = _beta_divergence(
            X, W, H, self._beta_loss, square_root=True, wX=None,
            W_weights=W_weights, H_weights=H_weights
        )

        self.n_components_ = H.shape[0]
        self.components_ = H
        self.n_iter_ = n_iter

        return W
    
    def transform(self, X, W_weights=None, H_weights=None):
        check_is_fitted(self)
        X = self._validate_data(
            X, accept_sparse=("csr", "csc"), dtype=[np.float64, np.float32], reset=False
        )

        with config_context(assume_finite=True):
            W, *_ = self._fit_transform(X, H=self.components_, update_H=False, W_weights=W_weights, H_weights=H_weights)

        return W

    def _fit_transform(self, X, y=None, W=None, H=None, update_H=True, W_weights=None, H_weights=None):
        if W_weights is None and H_weights is None:
            return super()._fit_transform(X, y, W, H, update_H=update_H)
        check_non_negative(X, "NMF (input X)")

        print(f'Weights stats: Median: {np.median(W_weights)}, Sum: {np.sum(W_weights)}, Max: {np.max(W_weights)}')
        
        # check parameters
        self._check_params(X)

        if X.min() == 0 and self._beta_loss <= 0:
            raise ValueError(
                "When beta_loss <= 0 and X contains zeros, "
                "the solver may diverge. Please add small values "
                "to X, or use a positive beta_loss."
            )

        # initialize or check W and H
        W, H = self._check_w_h(X, W, H, update_H)
        # scale the regularization terms
        l1_reg_W, l1_reg_H, l2_reg_W, l2_reg_H = self._compute_regularization(X)

        if self.solver == "cd":
            raise NotImplementedError
        elif self.solver == "mu":
            W, H, n_iter, *_ = _fit_multiplicative_update(
                X,
                W,
                H,
                self._beta_loss,
                self.max_iter,
                self.tol,
                l1_reg_W,
                l1_reg_H,
                l2_reg_W,
                l2_reg_H,
                update_H,
                self.verbose,
                W_weights=W_weights,
                H_weights=H_weights

            )
        else:
            raise ValueError("Invalid solver parameter '%s'." % self.solver)

        if n_iter == self.max_iter and self.tol > 0:
            warnings.warn(
                "Maximum number of iterations %d reached. Increase "
                "it to improve convergence."
                % self.max_iter,
                ConvergenceWarning,
            )

        return W, H, n_iter


def _fit_multiplicative_update(
    X,
    W,
    H,
    beta_loss="frobenius",
    max_iter=200,
    tol=1e-4,
    l1_reg_W=0,
    l1_reg_H=0,
    l2_reg_W=0,
    l2_reg_H=0,
    update_H=True,
    verbose=0,
    W_weights=None,
    H_weights=None
):
    start_time = time.time()

    beta_loss = _beta_loss_to_float(beta_loss)

    # gamma for Maximization-Minimization (MM) algorithm [Fevotte 2011]
    if beta_loss < 1:
        gamma = 1.0 / (2.0 - beta_loss)
    elif beta_loss > 2:
        gamma = 1.0 / (beta_loss - 1.0)
    else:
        gamma = 1.0

    if W_weights is None:
        wX = X
    else:
        wX = get_wX(X, W_weights, H_weights)

    # used for the convergence criterion
    error_at_init = _beta_divergence(
        X, W, H, beta_loss, square_root=True, wX=wX,
        W_weights=W_weights, H_weights=H_weights
    )
    print('Error at init', error_at_init, flush=True)
    previous_error = error_at_init

    H_sum, HHt, XHt = None, None, None

    for n_iter in range(1, max_iter + 1):
        # update W
        # H_sum, HHt and XHt are saved and reused if not update_H
        W, H_sum, HHt, XHt = _multiplicative_update_w(
            wX,
            W,
            H,
            beta_loss=beta_loss,
            l1_reg_W=l1_reg_W,
            l2_reg_W=l2_reg_W,
            gamma=gamma,
            H_sum=H_sum,
            HHt=HHt,
            XHt=XHt,
            update_H=update_H,
            W_weights=W_weights,
            H_weights=H_weights
        )

        # necessary for stability with beta_loss < 1
        if beta_loss < 1:
            W[W < EPSILON] = 0.0

        # update H (only at fit or fit_transform)
        if update_H:
            H = _multiplicative_update_h(
                wX,
                W,
                H,
                beta_loss=beta_loss,
                l1_reg_H=l1_reg_H,
                l2_reg_H=l2_reg_H,
                gamma=gamma,
                W_weights=W_weights,
                H_weights=H_weights
            )

            # These values will be recomputed since H changed
            H_sum, HHt, XHt = None, None, None

            # necessary for stability with beta_loss < 1
            if beta_loss <= 1:
                H[H < EPSILON] = 0.0

        # test convergence criterion every 10 iterations
        if tol > 0 and n_iter % 10 == 0:
            error = _beta_divergence(
                X, W, H, beta_loss,
                square_root=True,
                wX=wX,
                W_weights=W_weights,
                H_weights=H_weights
            )

            if verbose:
                iter_time = time.time()
                print(
                    "Epoch %02d reached after %.3f seconds, error: %f"
                    % (n_iter, iter_time - start_time, error),
                    flush=True
                )

            if (previous_error - error) / error_at_init < tol:
                break
            previous_error = error

    # do not print if we have already printed in the convergence test
    if verbose and (tol == 0 or n_iter % 10 != 0):
        end_time = time.time()
        print(
            "Epoch %02d reached after %.3f seconds." % (n_iter, end_time - start_time),
            flush=True
        )

    return W, H, n_iter


def _multiplicative_update_w(
    X,
    W,
    H,
    beta_loss,
    l1_reg_W,
    l2_reg_W,
    gamma,
    H_sum=None,
    HHt=None,
    XHt=None,
    update_H=True,
    W_weights=None,
    H_weights=None
):
    """Update W in Multiplicative Update NMF."""
    if beta_loss == 2:
        # Numerator
        if XHt is None:
            XHt = safe_sparse_dot(X, H.T)
        if update_H:
            # avoid a copy of XHt, which will be re-computed (update_H=True)
            numerator = XHt
        else:
            # preserve the XHt, which is not re-computed (update_H=False)
            numerator = XHt.copy()

        # Denominator
        if HHt is None:
            if H_weights is not None:
                HHt = np.dot(H_weights * H, H.T)
            else:
                HHt = np.dot(H, H.T)

        if W_weights is None:
            denominator = np.dot(W, HHt)
        else:
            denominator = np.dot(W_weights * W, HHt)
    else:
        raise NotImplementedError

    # Add L1 and L2 regularization
    if l1_reg_W > 0:
        denominator += l1_reg_W
    if l2_reg_W > 0:
        denominator = denominator + l2_reg_W * W
    denominator[denominator == 0] = EPSILON

    numerator /= denominator
    delta_W = numerator

    # gamma is in ]0, 1]
    if gamma != 1:
        delta_W **= gamma

    W *= delta_W

    return W, H_sum, HHt, XHt


def _multiplicative_update_h(
    X, W, H, beta_loss, l1_reg_H, l2_reg_H, gamma, A=None, B=None, rho=None,
    W_weights=None, H_weights=None
):
    """update H in Multiplicative Update NMF."""
    if beta_loss == 2:
        numerator = safe_sparse_dot(W.T, X)
        if W_weights is None:
            denominator = np.linalg.multi_dot([W.T, W, H])
        else:
            denominator = np.linalg.multi_dot([W.T, W_weights * W, H * H_weights])

    else:
        raise NotImplementedError

    # Add L1 and L2 regularization
    if l1_reg_H > 0:
        denominator += l1_reg_H
    if l2_reg_H > 0:
        denominator = denominator + l2_reg_H * H
    denominator[denominator == 0] = EPSILON

    if A is not None and B is not None:
        # Updates for the online nmf
        if gamma != 1:
            H **= 1 / gamma
        numerator *= H
        A *= rho
        B *= rho
        A += numerator
        B += denominator
        H = A / B

        if gamma != 1:
            H **= gamma
    else:
        delta_H = numerator
        delta_H /= denominator
        if gamma != 1:
            delta_H **= gamma
        H *= delta_H

    return H


def get_wX(X, W_weights, H_weights):
    if sp.issparse(X):
        wX = X.multiply(H_weights.squeeze()).T.multiply(W_weights.squeeze()).T.tocsr()
    else:
        wX = X * W_weights * H_weights
    return wX


def _beta_divergence(X, W, H, beta, square_root=False, wX=None, W_weights=None, H_weights=None):
    beta = _beta_loss_to_float(beta)

    # The method can be called with scalars
    if not sp.issparse(X):
        X = np.atleast_2d(X)
    W = np.atleast_2d(W)
    H = np.atleast_2d(H)

    # Frobenius norm
    if beta == 2:
        # Avoid the creation of the dense np.dot(W, H) if X is sparse.
        if sp.issparse(X):
            if W_weights is None:
                norm_X = np.dot(X.data, X.data)
                norm_WH = trace_dot(np.linalg.multi_dot([W.T, W, H]), H)
                cross_prod = trace_dot((X @ H.T), W)
            else:
                if wX is None:
                    wX = get_wX(X, W_weights, H_weights)
                sqrt_wW = W * np.sqrt(W_weights)
                sqrt_wH = H * np.sqrt(H_weights)

                norm_X = np.dot(wX.data, X.data)
                norm_WH = trace_dot(np.linalg.multi_dot([sqrt_wW.T, sqrt_wW, sqrt_wH]), sqrt_wH)
                cross_prod = trace_dot((wX @ H.T), W)

            res = (norm_X + norm_WH - 2.0 * cross_prod) / 2.0

        else:
            if W_weights is None:
                res = squared_norm(X - np.dot(W, H)) / 2.0
            else:
                res = squared_norm(
                    (X - np.dot(W, H)) * np.sqrt(W_weights) * np.sqrt(H_weights)
                ) / 2.0

        if square_root:
            if res < 0:
                print(res, norm_X, norm_WH, cross_prod)
            res = max(res, 0)
            return np.sqrt(res * 2)
        else:
            return res
    else:
        raise NotImplementedError
