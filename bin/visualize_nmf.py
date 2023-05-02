import os
import sys
import time
import numpy as np
import pandas as pd
import gzip
import matplotlib.pyplot as plt
import umap
import glob


import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.font_manager as font_manager
import seaborn as sns
import random
import math


from sklearn.decomposition._nmf import _BaseNMF, _beta_loss_to_float, _initialize_nmf, _fit_coordinate_descent
from abc import ABC
from numbers import Integral, Real
import numpy as np
import scipy.sparse as sp
import time
import itertools
import warnings
from math import sqrt
from scipy import linalg


#from ._cdnmf_fast import _update_cdnmf_fast
from sklearn._config import config_context
from sklearn.base import BaseEstimator, TransformerMixin, ClassNamePrefixFeaturesOutMixin
from sklearn.exceptions import ConvergenceWarning
from sklearn.utils import check_random_state, check_array, gen_batches
from sklearn.utils.extmath import randomized_svd, safe_sparse_dot, squared_norm
from sklearn.utils.validation import (
    check_is_fitted,
    check_non_negative,
)
from sklearn.utils._param_validation import (
    Interval,
    StrOptions,
    validate_params,
)

sys.path.append('/home/sboytsov/Vocabulary')
import OONMF
import OONMFhelpers
from OONMFhelpers import get_barsortorder

from common import get_matrix_path

