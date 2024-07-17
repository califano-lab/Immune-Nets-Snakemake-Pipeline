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
barcodes = pd.read_csv(metadata_file, header=None, sep='\t')
features = pd.read_csv(cluster_file, header=None, sep='\t')

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
sc.pl.umap(gene_expr_signature, color='leiden')



