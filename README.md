# Project Title: Analyzing scRNA-seq Data with ARACNe3 and PyViper

This repository contains a pipeline for analyzing single-cell RNA-seq (scRNA-seq) data. The analysis includes preprocessing, quality control, normalization, gene regulatory network inference using ARACNe3, and regulon activity analysis using PyViper.

## Pipeline Overview
The pipeline is managed using a Snakefile, which defines the rules for each step in the analysis. Below is the content of the Snakefile, which outlines how the data is processed and how each step is connected:

```
rule all:
    input:
        heatmap="figures/heatmap.png",
        umap_plot="figures/umap.png"

rule load_and_preprocess:
    input:
        gene_expr="Tutorial_1_gExpr_fibroblast_5802.tsv",
        network="fibroblast-net.tsv"

    output:
        processed_expr="results/processed_expr.h5ad",
        processed_net="results/processed_net.pkl"
    script:
        "scripts/load_and_preprocess.py"

rule viper_and_pca_analysis:
    input:
        processed_expr="results/processed_expr.h5ad",
        processed_net="results/processed_net.pkl"
    output:
        prot_act_pca="results/prot_act_pca.h5ad"
    script:
        "scripts/viper_and_pca_analysis.py"

rule clustering_and_umap:
    input:
        prot_act_pca="results/prot_act_pca.h5ad"
    output:
        umap_data="results/umap_data.h5ad",
        umap_plot="figures/umap.png"
    script:
        "scripts/clustering_and_umap.py"

rule integration_and_heatmap:
    input:
        umap_data="results/umap_data.h5ad"
    output:
        heatmap="figures/heatmap.png"
    script:
        "scripts/integration_and_heatmap.py"

```

## Description of Snakefile Rules

1. all:

Specifies the final outputs of the entire workflow, i.e., the UMAP plot (umap_plot) and the heatmap (heatmap).
These will be generated as a result of running the entire pipeline.

2. load_and_preprocess:

Input: Takes the gene expression data (Tutorial_1_gExpr_fibroblast_5802.tsv) and the network data (fibroblast-net.tsv).

Output: Produces preprocessed gene expression data (processed_expr.h5ad) and a processed network (processed_net.pkl).

Script: Runs the script load_and_preprocess.py to carry out this task.

3. viper_and_pca_analysis:

Input: Takes the processed expression data and network from the previous step.

Output: Outputs the protein activity PCA data (prot_act_pca.h5ad).

Script: Executes viper_and_pca_analysis.py to perform VIPER analysis and PCA.

4. clustering_and_umap:

Input: Uses the PCA data from the previous step.

Output: Outputs the UMAP plot data (umap_data.h5ad) and the UMAP plot image (umap_plot.png).

Script: Runs the clustering_and_umap.py script for UMAP and clustering.

5. integration_and_heatmap:

Input: Takes the UMAP data from the previous step.

Output: Produces the final heatmap (heatmap.png).

Script: Executes integration_and_heatmap.py to generate the heatmap.


## Introduction

This pipeline is designed for the analysis of scRNA-seq data to explore gene regulatory networks and protein activity levels.

## Installation

To install the necessary dependencies, run:

```bash
pip install scanpy pyviper
```

##  1. Preprocessing Data

This section handles the preprocessing of raw gene expression data and metadata to create an AnnData object that will be used for further analysis.

```import pandas as pd
import scanpy as sc
import pyviper
import pickle

def load_and_preprocess(gene_expr_path, network_path, output_expr, output_net):
    gene_expr_signature = pd.read_csv(gene_expr_path, sep="\t", index_col=0)
    gene_expr_signature = sc.AnnData(gene_expr_signature)
    
    network = pd.read_csv(network_path, delimiter="\t")
    network_interactome = pyviper.Interactome('immune', network)
    network_interactome.filter_targets(gene_expr_signature.var_names)
    
    # Save processed files
    gene_expr_signature.write_h5ad(output_expr)
    with open(output_net, 'wb') as f:
        pickle.dump(network_interactome, f)

if __name__ == "__main__":
    load_and_preprocess(snakemake.input.gene_expr, snakemake.input.network, snakemake.output.processed_expr, snakemake.output.processed_net)
```

##  2. Viper and pca analysis
This step performs VIPER analysis to estimate protein activity from gene expression data and then applies PCA for dimensionality reduction.
```
import scanpy as sc
import pyviper
import pickle

def run_viper_and_pca(processed_expr_path, network_interactome_path):
    # Load data
    gene_expr_signature = sc.read_h5ad(processed_expr_path)
    with open(network_interactome_path, 'rb') as f:
        network_interactome = pickle.load(f)
    
    # Perform VIPER analysis and PCA
    ProtAct_NaRnEA = pyviper.viper(gex_data=gene_expr_signature, interactome=network_interactome, enrichment="narnea", eset_filter=False, njobs=1, verbose=False)
    pyviper.tl.pca(ProtAct_NaRnEA, layer="pes", zero_center=True, svd_solver='arpack', random_state=0)
    return ProtAct_NaRnEA

if __name__ == "__main__":
    ProtAct_NaRnEA = run_viper_and_pca(snakemake.input.processed_expr, snakemake.input.processed_net)
    ProtAct_NaRnEA.write_h5ad(snakemake.output.prot_act_pca)
```

##  3. UMAP and Clustering
This step performs UMAP for dimensionality reduction and clustering analysis using the Leiden algorithm.
```
import scanpy as sc
import matplotlib.pyplot as plt
import os

def perform_clustering_and_umap(input_file, output_file, output_figure):
    data = sc.read_h5ad(input_file)
    
    # Perform neighborhood analysis, clustering, and UMAP
    sc.pp.neighbors(data, n_neighbors=20, n_pcs=50)
    sc.tl.leiden(data, resolution=0.1)
    sc.tl.umap(data)

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_figure), exist_ok=True)

    # Save UMAP plot
    sc.pl.umap(data, color='leiden', show=False, save=False)  # 先生成图，不直接保存
    plt.savefig(output_figure)  # 使用matplotlib保存，确保路径正确
    plt.close()

    # Save processed data
    data.write_h5ad(output_file)

if __name__ == "__main__":
    perform_clustering_and_umap(
        snakemake.input.prot_act_pca, 
        snakemake.output.umap_data, 
        snakemake.output.umap_plot
    )

```
##  4. Heatmap Generation
generates heatmaps to visualize the top activated proteins across clusters.
```
import pandas as pd
import scanpy as sc
import pyviper
import pickle

def load_and_preprocess(gene_expr_path, network_path, output_expr, output_net):
    gene_expr_signature = pd.read_csv(gene_expr_path, sep="\t", index_col=0)
    gene_expr_signature = sc.AnnData(gene_expr_signature)
    
    network = pd.read_csv(network_path, delimiter="\t")
    network_interactome = pyviper.Interactome('immune', network)
    network_interactome.filter_targets(gene_expr_signature.var_names)
    
    # Save processed files
    gene_expr_signature.write_h5ad(output_expr)
    with open(output_net, 'wb') as f:
        pickle.dump(network_interactome, f)

if __name__ == "__main__":
    load_and_preprocess(snakemake.input.gene_expr, snakemake.input.network, snakemake.output.processed_expr, snakemake.output.processed_net)
```
## Usage
To run the entire workflow, navigate to the directory containing the Snakefile and execute the following command:
bash code:
```
snakemake --cores 1 --snakefile Snakefile2
```


![umap](https://github.com/user-attachments/assets/599f3d12-6770-4a7d-81fa-22217a25488e)
![heatmap](https://github.com/user-attachments/assets/d8ce1c10-7ab1-45cd-ba87-53205146e09f)


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




