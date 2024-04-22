import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import anndata as ad
import json
from tqdm import tqdm
import os

def clr_mtx(t):

    y = t.copy()
    for i in range(y.shape[1]):
        x = t[:,i]
        positive_values = x[x > 0]
        log1p_positive = np.log1p(positive_values)

        # Calculate the sum of log1p values, ignoring NaNs
        sum_log1p = np.nansum(log1p_positive)

        # Calculate the geometric mean (in log space)
        geometric_mean_log_space = sum_log1p / len(x)
        x = np.where(x <= 0, np.finfo(float).tiny, x)

        # Apply log1p to the CLR transformation
        clr_transformed = np.log1p(x / np.exp(geometric_mean_log_space))
        y[:,i] = clr_transformed

    return y

if __name__ == '__main__':
    output_dir = 'output'
    adata = ad.read_h5ad(os.path.join(output_dir, 'raw.h5ad'))

    adata.layers['counts'] = adata.X.copy()
    adata.X = clr_mtx(adata.X)
    adata.write(os.path.join(output_dir, 'normalized_clr.h5ad'))