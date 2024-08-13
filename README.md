# Project Title: Analyzing scRNA-seq Data with ARACNe3 and PyViper

This repository contains a pipeline for analyzing single-cell RNA-seq (scRNA-seq) data. The analysis includes preprocessing, quality control, normalization, gene regulatory network inference using ARACNe3, and regulon activity analysis using PyViper.

## Table of Contents

1. Preprocessing Data
2. Quality Control
3. Normalization
4. PCA Analysis
5. UMAP Clustering
6. Differential Expression Analysis
7. ARACNe3 for Network Inference
8. Protein Activity Analysis with pyviper

## Introduction

This pipeline is designed for the analysis of scRNA-seq data to explore gene regulatory networks and protein activity levels.

## Installation

To install the necessary dependencies, run:

```bash
pip install scanpy pyviper
```

##  1. Preprocessing Data

This section handles the preprocessing of raw gene expression data and metadata to create an AnnData object that will be used for further analysis.

```import scanpy as sc
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
```

##  2. Normalizing Data
```After preprocessing, the data needs to be normalized. This section covers the steps to normalize and scale the data, making it ready for PCA and UMAP analysis.

import scanpy as sc

data_path = "/Users/lzy/Desktop/"
adata = sc.read(data_path + "preprocessed_data.h5ad")

# Normalization and transformation
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var['highly_variable']]
sc.pp.regress_out(adata, ['total_counts', 'percent_mito'])
sc.pp.scale(adata, max_value=10)

# Save normalized data
adata.write(data_path + "normalized_data.h5ad")
```

##  3. PCA Analysis
```
import scanpy as sc

data_path = "/Users/lzy/Desktop/"
adata = sc.read(data_path + "normalized_data.h5ad")

# PCA analysis
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color='Cell_Type')

# Save PCA data
adata.write(data_path + "pca_data.h5ad")
```
##  4. UMAP Clustering
Use UMAP for clustering analysis and visualize the results.
```
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
```

## 5. Differential Expression Analysis
Identify differentially expressed genes across clusters.
```
import scanpy as sc

data_path = "/Users/lzy/Desktop/"
adata = sc.read(data_path + "umap_data.h5ad")

# Differential expression analysis
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

# Save differential expression data
adata.write(data_path + "differential_expression_data.h5ad")
```

##  ARACNe-3 Analysis

1. Installing and load the human_regulators.rda file in R. You can achieve this by using the following R script:
```
# Install necessary packages if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Load the human_regulators.rda file
load("path_to_rda_file/human_regulators.rda")
```
2. Export the Four Regulators Files

```
# Save each regulators list as a CSV file
write.csv(human_regulators_tf, "human_regulators_tf.csv", row.names = FALSE)
write.csv(human_regulators_cotf, "human_regulators_cotf.csv", row.names = FALSE)
write.csv(human_regulators_sig, "human_regulators_sig.csv", row.names = FALSE)
write.csv(human_regulators_surf, "human_regulators_surf.csv", row.names = FALSE)
```

3. Use the following Python script to match the genes in the gene expression data with those in the regulators files:

```
import pandas as pd

# Load gene expression data
gene_expression = pd.read_csv('/path_to_your_file/gene_expression.csv', index_col=0)

# Convert gene names to upper case
genes_in_expression_data = gene_expression.index.str.upper().tolist()

# Load regulators files
tf = pd.read_csv("/path_to_your_file/human_regulators_tf.csv")
cotf = pd.read_csv("/path_to_your_file/human_regulators_cotf.csv")
sig = pd.read_csv("/path_to_your_file/human_regulators_sig.csv")
surf = pd.read_csv("/path_to_your_file/human_regulators_surf.csv")

# Convert regulator gene names to upper case
tf_genes = tf.iloc[:, 0].str.upper().tolist()
cotf_genes = cotf.iloc[:, 0].str.upper().tolist()
sig_genes = sig.iloc[:, 0].str.upper().tolist()
surf_genes = surf.iloc[:, 0].str.upper().tolist()

# Match regulators with gene expression data
matched_tf = [gene for gene in tf_genes if gene in genes_in_expression_data]
matched_cotf = [gene for gene in cotf_genes if gene in genes_in_expression_data]
matched_sig = [gene for gene in sig_genes if gene in genes_in_expression_data]
matched_surf = [gene for gene in surf_genes if gene in genes_in_expression_data]

# Combine all matched regulators and remove duplicates
all_matched_regulators = list(set(matched_tf + matched_cotf + matched_sig + matched_surf))

# Save the matched regulators to a CSV file
pd.DataFrame(all_matched_regulators, columns=["Regulators"]).to_csv("/path_to_your_file/matched_regulators.csv", index=False)
```
4. Generate Gene Regulatory Network Using ARACNe-3
use ARACNe-3 to generate the gene regulatory network by following these steps based on the five examples provided with ARACNe-3:

Example 1: Using Default Parameters
```./ARACNe3_app_release -e /path_to_your_file/gene_expression.csv -r /path_to_your_file/matched_regulators.csv -o /path_to_your_file/aracne_output/ --subsample 0.632 --threads 4```

Example 2: Specifying Edge P-Value Threshold
```./ARACNe3_app_release -e /path_to_your_file/gene_expression.csv -r /path_to_your_file/matched_regulators.csv -o /path_to_your_file/aracne_output/ --pvalue 1E-8 --subsample 0.632 --threads 4```

Example 3: Setting a Higher Subsampling Rate
```./ARACNe3_app_release -e /path_to_your_file/gene_expression.csv -r /path_to_your_file/matched_regulators.csv -o /path_to_your_file/aracne_output/ --subsample 0.8 --threads 4```

Example 4: Using Mutual Information Method
```./ARACNe3_app_release -e /path_to_your_file/gene_expression.csv -r /path_to_your_file/matched_regulators.csv -o /path_to_your_file/aracne_output/ --method mi --subsample 0.632 --threads 4```

Example 5: Changing the Number of Bootstrap Iterations
```./ARACNe3_app_release -e /path_to_your_file/gene_expression.csv -r /path_to_your_file/matched_regulators.csv -o /path_to_your_file/aracne_output/ --subsample 0.632 --threads 4 --bootstrap 100```

## Pyviper Analysis
Analyzing scRNA-seq data at the Protein Activity Level

1.Import modules
```import pyviper
import scanpy as sc
import anndata 
import pandas as pd
import numpy as np
import random```

```data_location = "/Users/lzy/Desktop/"
```

Step 1. Load a gene expression signature for single-cells
```gene_expr_path = data_location + "gene_expression_TPM"
gene_expr_signature = pd.read_csv(gene_expr_path, sep="\t")
gene_expr_signature = sc.AnnData(gene_expr_signature)
```
```gene_expr_signature```

Step 2. Load the gene regulatory network

```network_path = data_location + "TPM-net.tsv" # path to ARACNe network
    
network = pd.read_csv(network_path, delimiter="\t")```

```network_interactome = pyviper.Interactome('fibroblasts', network)
network_interactome.size()```

```network_interactome.net_table.head()```

Filter out targets in the Interactome
```network_interactome.filter_targets(gene_expr_signature.var_names)```

As an example, display the number of targets of two regulators, TSPAN6 and FGR.
```
n_TSPAN6 = len(network_interactome.get_reg('TSPAN6'))
n_FGR = len(network_interactome.get_reg('FGR'))

print("The number of targets of TSPAN6 and FGR are " + str(n_TSPAN6) + " and " + str(n_FGR) + " respectively.")```

Step 3. Convert the gene expression signature into a protein activity matrix using VIPER
```network_pruned = network_interactome.copy()
network_pruned.prune(max_targets=50,eliminate=True)```

Now all the regulators in the network have exactly 50 transcriptional targets. See e.g. TSPAN6 and FGR.
```n_TSPAN6 = len(network_pruned.get_reg('TSPAN6'))  # number of TSPAN6 targets in the network
n_FGR = len(network_pruned.get_reg('FGR'))  # number of FGR targets in the network

print("Number of TSPAN6 targets: " + str(n_TSPAN6) + "\nNumber of FGR targets: " + str(n_FGR))```

```ProtAct_aREA = pyviper.viper(gex_data=gene_expr_signature, # gene expression signature
                             interactome=network_pruned, # gene regulatory network
                             enrichment = "area",
                             output_as_anndata=False,
                             njobs=1,
                             verbose=False)```

ProtAct_aREA contains the activity of every regulatory protein described by its NES for each single cell.
```ProtAct_aREA```

Step 4. Analyze single-cells at the Protein Activity level

```pyviper.config.set_regulators_species_to_use(species="human")```




