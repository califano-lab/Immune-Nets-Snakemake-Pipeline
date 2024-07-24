import scanpy as sc
import numpy as np

data_path = "/Users/lzy/Desktop/"
adata = sc.read(data_path + "preprocessed_data.h5ad")

# Data preprocessing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.obs['percent_mito'] = np.sum(
    adata[:, adata.var_names.str.startswith('MT-')].X, axis=1) / np.sum(adata.X, axis=1)
adata.obs['total_counts'] = np.sum(adata.X, axis=1)
adata = adata[adata.obs['percent_mito'] < 0.05]

# Save quality controlled data
adata.write(data_path + "qc_data.h5ad")