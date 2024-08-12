import scanpy as sc

# Set the data path
data_path = "/Users/lzy/Desktop/"

# Load the UMAP data
adata = sc.read(data_path + "umap_data.h5ad")

# Perform differential expression analysis
sc.tl.rank_genes_groups(
    adata, 
    groupby='leiden',  # Group by clusters identified by the Leiden algorithm
    method='wilcoxon',  # Use Wilcoxon rank-sum test for differential expression
    key_added='rank_genes_leiden'  # Store the results in a specific key
)

# Plot the top differentially expressed genes for each cluster
sc.pl.rank_genes_groups(adata, key='rank_genes_leiden', n_genes=20, sharey=False)

# Save the differential expression results
adata.write(data_path + "differential_expression_data.h5ad")
