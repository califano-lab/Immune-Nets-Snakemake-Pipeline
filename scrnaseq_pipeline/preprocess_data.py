import scanpy as sc
import pandas as pd
import numpy as np

# Define file paths
data_path = "/Users/lzy/Desktop/"
matrix_file = data_path + "bm_scp_gex_matrix.csv"
meta_file = data_path + "bm_scp_meta_(1).txt"

# Load gene expression matrix
gene_expression = pd.read_csv(matrix_file, index_col=0)
gene_expression = gene_expression.T  # Transpose to make genes in columns

# Load metadata information
meta_data = pd.read_csv(meta_file, sep="\t", index_col=0)

# Check and align metadata and gene expression data
common_indices = gene_expression.index.intersection(meta_data.index)
gene_expression = gene_expression.loc[common_indices]
meta_data = meta_data.loc[common_indices]

# Create AnnData object
adata = sc.AnnData(gene_expression)
adata.obs = meta_data

# Save preprocessed data
adata.write(data_path + "preprocessed_data.h5ad")
