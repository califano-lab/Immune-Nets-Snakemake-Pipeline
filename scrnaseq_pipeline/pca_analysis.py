import scanpy as sc

data_path = "/Users/lzy/Desktop/"
adata = sc.read(data_path + "normalized_data.h5ad")

# PCA analysis
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color='Cell_Type')

# Save PCA data
adata.write(data_path + "pca_data.h5ad")