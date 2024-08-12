import scanpy as sc

# Set the data path
data_path = "/Users/lzy/Desktop/"

# Load the PCA data
adata = sc.read(data_path + "pca_data.h5ad")

# Perform clustering analysis and UMAP visualization
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)  # Compute the neighbor graph using PCA components
sc.tl.leiden(adata)  # Perform Leiden clustering
sc.tl.umap(adata)  # Perform UMAP dimensionality reduction
sc.pl.umap(adata, color=['Cell_Type', 'leiden'])  # Visualize UMAP colored by cell type and Leiden clusters

# Save the UMAP results
adata.write(data_path + "umap_data.h5ad")

