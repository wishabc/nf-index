import sys
import pandas as pd
import numpy as np
import os
from helpers import read_zarr_backed, convert_to_sparse_if_sufficently_sparse


def main(adata, meta, matrices):
    