import scanpy as sc

# Set the data path
data_path = "/Users/lzy/Desktop/"

# Load the normalized data
adata = sc.read(data_path + "normalized_data.h5ad")

# Perform PCA analysis
sc.tl.pca(adata, svd_solver='arpack')

# Visualize the PCA results
sc.pl.pca(adata, color='Cell_Type')

# Save the PCA results
adata.write(data_path + "pca_data.h5ad")
