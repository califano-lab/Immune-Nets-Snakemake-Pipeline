import pyviper
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import random
import warnings

warnings.filterwarnings("ignore", message=".*The 'nopython' keyword.*")  # for jit decorator issue with sc.pp.neighbors

# Set data paths
data_location = "/Users/lzy/Desktop/"
matrix_file = data_location + "expression_matrix.mtx.gz"
metadata_file = data_location + "PancrImmune1-metadata_modified.tsv"  # Use the modified metadata file
cluster_file = data_location + "PancrImmune1-cluster.tsv"

# Read gene expression signature
gene_expr_signature = sc.read_mtx(matrix_file).T  # Transpose matrix to have cells as rows and genes as columns
barcodes = pd.read_csv(metadata_file, header=None, sep='\t', low_memory=False)
features = pd.read_csv(cluster_file, header=None, sep='\t', low_memory=False)

# Debugging: Print lengths
print(f"Length of barcodes: {len(barcodes)}")
print(f"Length of gene_expr_signature.obs: {gene_expr_signature.obs.shape[0]}")
print(f"Length of features: {len(features)}")
print(f"Length of gene_expr_signature.var: {gene_expr_signature.var.shape[0]}")

# Print the first few rows of each to inspect
print(barcodes.head())
print(features.head())

# Ensure lengths match
if len(barcodes) != gene_expr_signature.obs.shape[0]:
    raise ValueError(f"Length of barcodes ({len(barcodes)}) does not match length of gene_expr_signature.obs ({gene_expr_signature.obs.shape[0]}).")

# Ensure that the number of genes in features matches the number of genes in the expression matrix
if len(features) != gene_expr_signature.shape[1]:
    raise ValueError(f"The number of genes in the features file ({len(features)}) does not match the number of genes in the expression matrix ({gene_expr_signature.shape[1]}).")

# Add barcodes and gene information
gene_expr_signature.obs['barcode'] = barcodes[0].values
gene_expr_signature.var['gene'] = features[1].values
gene_expr_signature.var_names = gene_expr_signature.var['gene']
gene_expr_signature.obs_names = gene_expr_signature.obs['barcode']

print(gene_expr_signature)

# Compute nearest neighbors
sc.pp.neighbors(gene_expr_signature, n_neighbors=10)

# Compute UMAP
sc.tl.umap(gene_expr_signature)

# Plot UMAP
sc.pl.umap(gene_expr_signature, color=['donor', 'sex', 'cell_type__ontology_label'])

# Cluster using Leiden algorithm
sc.tl.leiden(gene_expr_signature, resolution=0.1)

# Plot clustering result
sc.pl.umap(gene_expr_signature, color='leiden')
