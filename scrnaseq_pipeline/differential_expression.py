import scanpy as sc

data_path = "/Users/lzy/Desktop/"
adata = sc.read(data_path + "umap_data.h5ad")

# Differential expression analysis
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)

# Save differential expression data
adata.write(data_path + "final_data.h5ad")
