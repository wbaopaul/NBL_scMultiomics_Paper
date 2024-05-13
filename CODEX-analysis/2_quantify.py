import numpy as np
import scanpy as sc
import pandas as pd
import anndata as ad
from skimage import io
import os
from tqdm import tqdm
import numba

markers = np.loadtxt('', dtype=str)
img_path = ''
output_dir = 'output'
mask = io.imread(os.path.join(output_dir, 'cell_mesmer.tif'))


img_data = io.imread(img_path)
img_data = np.transpose(img_data, (1,2,0))

n_cells = np.max(mask)

# coords = np.array(np.ogrid[tuple(slice(dim) for dim in mask.shape)])
exp_rst = np.zeros((n_cells, img_data.shape[2]))
centers = np.zeros((n_cells, 2))
cell_ids = np.zeros((n_cells, 1))
for i in tqdm(range(n_cells)):
    exp_rst[i,:] = np.mean(img_data[mask==i+1,:], axis=0)
    centers[i,:] = np.mean(np.argwhere(mask==i+1), axis=0)
    cell_ids[i,:] = i+1

df = pd.DataFrame(exp_rst, columns=markers)
df['x'] = centers[:,0]
df['y'] = centers[:,1]
df['cell_id'] = cell_ids[:,0]
df.to_csv(os.path.join(output_dir, 'cell_gene.csv'), index=False)

exp_rst = df[markers].values
centers = df[['x', 'y']].values
adata = ad.AnnData(X=exp_rst, obs=pd.DataFrame(index=['cell_{}'.format(i) for i in range(n_cells)]), var=pd.DataFrame(index=markers))
adata.obsm['spatial'] = centers
adata.obs['cell_id'] = cell_ids[:,0]
print(adata)
adata.write(os.path.join(output_dir, 'raw.h5ad'))




