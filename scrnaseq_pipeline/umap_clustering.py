import scanpy as sc

data_path = "/Users/lzy/Desktop/"
adata = sc.read(data_path + "pca_data.h5ad")

# Clustering analysis and UMAP visualization
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color=['Cell_Type', 'leiden'])

# Save UMAP data
adata.write(data_path + "umap_data.h5ad")
