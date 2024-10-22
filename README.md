# Project Title: scRNA-seq Data Analysis Pipeline with PyVIPER

This repository contains a pipeline for analyzing single-cell RNA-seq (scRNA-seq) data. The analysis includes preprocessing, quality control, normalization, gene regulatory network inference using ARACNe3, and regulon activity analysis using PyViper.

## Overview
The goal of this pipeline is to process gene expression datasets and analyze them to understand regulatory networks and protein activity in single cells. This process involves multiple steps, including data preparation, ARACNe3 network generation, VIPER analysis, PCA, clustering, and visualization.

## Introduction of Snakemake

Snakemake is a workflow management system that allows you to create reproducible and scalable data analyses. Inspired by the Python programming language and the "Make" utility commonly used for build automation, Snakemake enables you to define complex workflows using simple rules that specify how input files are converted into output files.

## Pipeline Overview

```
# Read input and output folder from config or command-line arguments
input_folder = config.get("input_folder", "data")
output_folder = config.get("output_folder", "results")

# Detect all gene expression files in the input folder dynamically
import os
import glob

gene_expr_files = sorted(glob.glob(f"{input_folder}/*.tsv"))
datasets = [os.path.splitext(os.path.basename(f))[0] for f in gene_expr_files]

rule all:
    input:
        expand(f"{output_folder}/{{dataset}}_heatmap.png", dataset=datasets),
        expand(f"{output_folder}/{{dataset}}_umap.png", dataset=datasets),
        expand(f"{output_folder}/{{dataset}}_silhouette_plot.png", dataset=datasets),
        expand(f"{output_folder}/{{dataset}}_resolution_silhouette_plot.png", dataset=datasets),

# Step 2: Run ARACNe3 to generate networks
rule run_aracne3:
    input:
        expr_matrix=f"{input_folder}/{{dataset}}.tsv",
        regulators="combined_regulators.txt"
    output:
        directory(f"{output_folder}/{{dataset}}_consolidated-net_defaultid")
    shell:
        """
        /Users/lzy/Desktop/ARACNe3/build/src/app/ARACNe3_app_release \
        -e {input.expr_matrix} \
        -r {input.regulators} \
        -o {output}/{{wildcards.dataset}}_consolidated-net_defaultid \
        -x 10 --alpha 0.05 --threads 1
        """

# Step 3: Convert ARACNe3 output to formatted networks
rule format_network:
    input:
        aracne_output=f"{output_folder}/{{dataset}}_consolidated-net_defaultid/consolidated-net_defaultid.tsv"
    output:
        formatted_network=f"{output_folder}/{{dataset}}_consolidated-net_defaultid_formatted_network.tsv"
    script:
        "scripts/format_network.py"

# Step 4: Load and preprocess data for VIPER
rule load_and_preprocess:
    input:
        gene_expr=f"{input_folder}/{{dataset}}.tsv",
        formatted_network=f"{output_folder}/{{dataset}}_consolidated-net_defaultid_formatted_network.tsv"
    output:
        processed_expr=f"{output_folder}/{{dataset}}_processed_expr.h5ad",
        processed_net=f"{output_folder}/{{dataset}}_processed_net.pkl"
    script:
        "scripts/load_and_preprocess.py"

# Step 5: VIPER analysis and PCA
rule viper_and_pca_analysis:
    input:
        processed_expr=f"{output_folder}/{{dataset}}_processed_expr.h5ad",
        processed_net=f"{output_folder}/{{dataset}}_processed_net.pkl"
    output:
        prot_act_pca=f"{output_folder}/{{dataset}}_prot_act_pca.h5ad"
    conda:
        "environment.yml"  # Specify the conda environment
    script:
        "scripts/viper_and_pca_analysis.py"

# Step 6: Clustering and UMAP analysis
rule clustering_and_umap:
    input:
        prot_act_pca=f"{output_folder}/{{dataset}}_prot_act_pca.h5ad"
    output:
        umap_data=f"{output_folder}/{{dataset}}_umap_data.h5ad",
        umap_plot=f"{output_folder}/{{dataset}}_umap.png"
    conda:
        "environment.yml"
    script:
        "scripts/clustering_and_umap.py"

# Step 7: Optimize resolution and generate silhouette plot for resolution vs score
rule resolution_optimization:
    input:
        processed_expr=f"{output_folder}/{{dataset}}_processed_expr.h5ad"
    output:
        silhouette_resolution_plot=f"{output_folder}/{{dataset}}_resolution_silhouette_plot.png",
        best_resolution=f"{output_folder}/{{dataset}}_best_resolution.txt"
    script:
        "scripts/resolution_optimization.py"

# Step 8: Generate Silhouette Plot
rule silhouette_plot:
    input:
        umap_data=f"{output_folder}/{{dataset}}_umap_data.h5ad",
        best_resolution=f"{output_folder}/{{dataset}}_best_resolution.txt"
    output:
        silhouette_plot=f"{output_folder}/{{dataset}}_silhouette_plot.png"
    script:
        "scripts/silhouette_plot.py"

# Step 9: Generate Heatmap
rule integration_and_heatmap:
    input:
        umap_data=f"{output_folder}/{{dataset}}_umap_data.h5ad"
    output:
        heatmap=f"{output_folder}/{{dataset}}_heatmap.png"
    script:
        "scripts/integration_and_heatmap.py"
```


## Description of Snakefile Rules

1. all:

Purpose: This is a meta-rule that defines the final outputs of the entire workflow. It ensures that all required steps have been completed and all outputs have been generated successfully.

Inputs:
· Heatmap (results/{dataset}_heatmap.png)
· UMAP plot (results/{dataset}_umap.png)
· Silhouette plot (results/{dataset}_silhouette_plot.png)
· Resolution-based silhouette plot (results/{dataset}_resolution_silhouette_plot.png)

Functionality: Acts as the endpoint to verify that all datasets have been processed through the pipeline.


2.run_aracne3:

Purpose: Executes ARACNe3, a tool for inferring gene regulatory networks, using gene expression data and a combined regulators file.

Inputs:

· Gene expression file (```data/{dataset}.tsv```).
· Combined regulators file (```combined_regulators.txt```).

Outputs:

A directory containing the ARACNe3 consolidated network output (```results/{dataset}_consolidated-net_defaultid```).

Shell Command: Runs the ARACNe3 executable with specified parameters to generate the network.

3.format_network:

Purpose: Converts the output from ARACNe3 into a format that is compatible with PyViper for subsequent analysis.

Inputs:
· ARACNe3 output file (```results/{dataset}_consolidated-net_defaultid/consolidated-net_defaultid.tsv```).

Outputs:
A formatted network file (```results/{dataset}_consolidated-net_defaultid_formatted_network.tsv```).

Script: Uses ```scripts/format_network.py``` to format the network data into a usable structure for PyViper.

4.load_and_preprocess:

Purpose: Loads gene expression data and formatted network data, and performs preprocessing necessary for downstream VIPER analysis.

Inputs:

Gene expression file (```data/{dataset}.tsv```).
Formatted network file (```results/{dataset}_consolidated-net_defaultid_formatted_network.tsv```).

Outputs:
· Preprocessed gene expression data (```results/{dataset}_processed_expr.h5ad```).
· Processed network file (```results/{dataset}_processed_net.pkl```).

Script: Executes ```scripts/load_and_preprocess.py``` to preprocess the data.

5. viper_and_pca_analysis:

Purpose: Conducts VIPER analysis to infer protein activity from gene expression data and applies PCA for dimensionality reduction.

Inputs:
Preprocessed gene expression data (```results/{dataset}_processed_expr.h5ad```).
Processed network file (```results/{dataset}_processed_net.pkl```).

Outputs:
Protein activity PCA data (```results/{dataset}_prot_act_pca.h5ad```).
Script: Runs ```scripts/viper_and_pca_analysis.py``` to perform VIPER and PCA analysis..

6. clustering_and_umap:

Purpose: Performs clustering analysis and generates UMAP visualizations to visualize the data in two-dimensional space.

Inputs:

Protein activity PCA data (```results/{dataset}_prot_act_pca.h5ad```)

Outputs:

UMAP data (```results/{dataset}_umap_data.h5ad```).
UMAP plot image (```results/{dataset}_umap.png```).

Script: Executes ```clustering_and_umap.py``` to perform clustering and generate UMAP plots.

7. resolution_optimization:

Purpose: Optimizes the clustering resolution by calculating silhouette scores across different resolution values and identifies the best resolution.

Inputs:

Processed expression data (results/{dataset}_processed_expr.h5ad)

Outputs:

Silhouette resolution plot (results/{dataset}_resolution_silhouette_plot.png)
Best resolution value (results/{dataset}_best_resolution.txt)

Script: Runs a Python script (scripts/resolution_optimization.py) to compute silhouette scores and determine the optimal clustering resolution.

8. silhouette_plot:

Purpose: Generates a silhouette plot based on the best resolution from the resolution_optimization step, providing a measure of cluster quality.

Inputs:

UMAP data (results/{dataset}_umap_data.h5ad)
Best resolution value (results/{dataset}_best_resolution.txt)

Outputs:

Silhouette plot (results/{dataset}_silhouette_plot.png)

Script: Runs a Python script (scripts/silhouette_plot.py) to create a silhouette plot based on the optimized clustering resolution.

9. integration_and_heatmap:

Purpose: Generates a heatmap based on the UMAP clustering and processed expression data to visualize gene expression across different cell types or clusters.

Inputs:
UMAP data (```results/{dataset}_umap_data.h5ad```).

Outputs:

Final heatmap image (```results/{dataset}_heatmap.png```).

Script: Runs ```scripts/integration_and_heatmap.py``` to create a heatmap visualization.


## Scripts Part



## 1. format_network
Shell code:
```
python python scripts/generate_h5ad.py

        /Users/lzy/Desktop/ARACNe3/build/src/app/ARACNe3_app_release \
        -e {input.expr_matrix} \
        -r {input.regulators} \
        -o results/{wildcards.dataset}_consolidated-net_defaultid \
        -x 10 --alpha 0.05 --threads 1

```
Python code:
```

import os
import pyviper
import pandas as pd

# Define the input and output root directories
input_root_dir = "results"
output_root_dir = "results"

# Retrieve all subdirectories, assuming all subdirectory names start with "gene_expression_sub_"
subdirs = sorted([d for d in os.listdir(input_root_dir) if d.startswith("gene_expression_sub_")])

for subdir in subdirs:
    input_file = os.path.join(input_root_dir, subdir, "consolidated-net_defaultid.tsv")
    output_file = os.path.join(output_root_dir, f"{subdir}_formatted_network.tsv")

    if os.path.exists(input_file):
        try:
            # Convert the ARACNe3 output to a regulon DataFrame
            regulon_df = pyviper.pp.aracne3_to_regulon(
                net_file=input_file,
                MI_thres=0,                   # Threshold for Mutual Information
                regul_size=50,                # Maximum number of targets per regulator
                normalize_MI_per_regulon=True # Normalize MI scores within each regulon
            )

            # Save the formatted network to a file
            regulon_df.to_csv(output_file, sep='\t', index=False)
            print(f"Formatted network for {subdir} saved to {output_file}.")
        except Exception as e:
            print(f"Failed to process {input_file} due to error: {e}")
    else:
        print(f"Input file {input_file} does not exist.")

```

## 2. load_and_preprocess


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
    load_and_preprocess(snakemake.input.gene_expr, snakemake.input.formatted_network, snakemake.output.processed_expr, snakemake.output.processed_net)
```

##  3. Viper and pca analysis

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

##  4. UMAP and Clustering

```
import scanpy as sc
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score
import os

def perform_clustering_and_umap(input_file, output_file, output_figure):
    # Load the input data
    data = sc.read_h5ad(input_file)
    
    # Perform neighborhood analysis, clustering, and UMAP
    sc.pp.neighbors(data, n_neighbors=15, n_pcs=50)
    sc.tl.leiden(data, resolution=0.5)
    sc.tl.umap(data)

    # Calculate silhouette score
    if 'leiden' in data.obs:
        silhouette_avg = silhouette_score(data.obsm['X_umap'], data.obs['leiden'])
        print(f"Silhouette Score: {silhouette_avg}")
    else:
        print("Leiden clustering not found in data.obs")

    # Save the UMAP data
    data.write_h5ad(output_file)
    
    # Ensure the output directory exists
    os.makedirs(os.path.dirname(output_figure), exist_ok=True)
    
    # Save the UMAP plot
    sc.pl.umap(data, color='leiden', show=False)  # Plot without saving
    
    # Save the plot using matplotlib
    plt.savefig(output_figure)  # Explicitly save the figure to the path
    plt.close()

if __name__ == "__main__":
    perform_clustering_and_umap(
        snakemake.input[0],  # Input data file
        snakemake.output[0],  # Output UMAP data
        snakemake.output[1]   # Output UMAP figure
    )

```



##  5. Optimize resolution and generate silhouette plot for resolution vs score

```
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score

def optimize_resolution(adata, resolutions, n_neighbors=15):
    scores = []
    for res in resolutions:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors)
        sc.tl.leiden(adata, resolution=res)
        
        if len(adata.obs['leiden'].unique()) > 1:
            score = silhouette_score(adata.obsm['X_pca'], adata.obs['leiden'])
            scores.append(score)
        else:
            scores.append(-1)  # Append -1 if only one cluster is found

    return scores

# Load the processed expression data
adata = sc.read_h5ad(snakemake.input.processed_expr)

# Define a range of resolutions to test
resolutions = np.linspace(0.1, 1.0, 10)

# Optimize silhouette score across resolutions
silhouette_scores = optimize_resolution(adata, resolutions)

# Find the best resolution based on the highest silhouette score
best_resolution = resolutions[np.argmax(silhouette_scores)]

# Save the best resolution to a file
with open(snakemake.output.best_resolution, 'w') as f:
    f.write(f"Best resolution: {best_resolution}\n")

# Plot the silhouette scores vs. resolution
plt.errorbar(resolutions, silhouette_scores, yerr=None, fmt='-o')
plt.xlabel('Resolution')
plt.ylabel('Silhouette Score')
plt.title(f"Best resolution = {best_resolution:.2f}")
plt.savefig(snakemake.output.silhouette_resolution_plot)
plt.close()
```

##  6. Generate Silhouette Plot

```
import scanpy as sc
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_samples, silhouette_score
import numpy as np
import os

def generate_silhouette_plot(input_file, output_file):
    # Load the input data
    data = sc.read_h5ad(input_file)

    # Ensure 'leiden' clusters exist
    if 'leiden' not in data.obs:
        print("Leiden clustering not found, performing clustering...")
        sc.pp.neighbors(data, n_neighbors=15, n_pcs=50)
        sc.tl.leiden(data, resolution=0.5)

    # Calculate silhouette scores for each cell
    if 'leiden' in data.obs:
        labels = data.obs['leiden'].astype(int)
        X = data.obsm['X_umap']  # UMAP embedding

        # Compute silhouette scores for each cell
        silhouette_vals = silhouette_samples(X, labels)

        # Calculate the average silhouette score
        silhouette_avg = np.mean(silhouette_vals)
        print(f"Average Silhouette Score: {silhouette_avg}")

        # Create silhouette plot for each cluster
        fig, ax = plt.subplots(figsize=(10, 6))

        y_lower, y_upper = 0, 0
        unique_labels = np.unique(labels)
        for i, label in enumerate(unique_labels):
            # Aggregate silhouette scores for the cluster
            ith_silhouette_vals = silhouette_vals[labels == label]
            ith_silhouette_vals.sort()

            size_cluster_i = ith_silhouette_vals.shape[0]
            y_upper = y_lower + size_cluster_i

            color = plt.cm.nipy_spectral(float(i) / len(unique_labels))
            ax.fill_betweenx(np.arange(y_lower, y_upper), 0, ith_silhouette_vals,
                             facecolor=color, edgecolor=color, alpha=0.7)

            ax.text(-0.05, y_lower + 0.5 * size_cluster_i, str(label))

            y_lower = y_upper  # Update y_lower for the next cluster

        ax.set_title(f"Silhouette Plot (Avg: {silhouette_avg:.3f})")
        ax.set_xlabel("Silhouette Coefficient")
        ax.set_ylabel("Cluster")

        # Draw vertical line for average silhouette score
        ax.axvline(x=silhouette_avg, color="red", linestyle="--")

        # Save the silhouette plot
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        plt.savefig(output_file)
        plt.close()

    else:
        print("Leiden clustering not available in the dataset.")

if __name__ == "__main__":
    # Replace these with Snakemake input/output or other parameters
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]
    
    generate_silhouette_plot(input_file, output_file)
```

##  7. Heatmap Generation

```
import scanpy as sc
import pyviper
import matplotlib.pyplot as plt
import numpy as np

def generate_heatmap(data_path, output_path, N=50):
    # Ensure all previous plot windows are closed
    plt.close('all')
    
    data = sc.read_h5ad(data_path)
    
    # Process data to ensure there are no invalid values
    data.X = np.where(data.X <= 0, np.min(data.X[data.X > 0]), data.X)
    
    try:
        # Select the top N most active proteins
        protein_set = data.var_names[:N]
        
        # Create Heatmap
        pyviper.pl.heatmap(data, var_names=protein_set, groupby="leiden", vcenter=0, cmap="RdBu_r", swap_axes=True, show=False)
        
        # Save plot
        plt.savefig(output_path, bbox_inches="tight")
        plt.close()
        
    except Exception as e:
        print(f"Failed to generate heatmap due to: {e}")

if __name__ == "__main__":
    generate_heatmap(snakemake.input[0], snakemake.output[0])

```
## Usage
To run the entire workflow, navigate to the directory containing the Snakefile and execute the following command:
bash code:
```
snakemake --cores 1 --snakefile Snakefile2

```
Example Visulization:

![gene_expression_sub_1_heatmap](https://github.com/user-attachments/assets/1f7c4850-1f9b-47ce-b0f8-b9c816c8804b)
![gene_expression_sub_1_resolution_silhouette_plot](https://github.com/user-attachments/assets/6e15029c-2b92-43e9-8e4f-2fc028cc198e)
![gene_expression_sub_1_silhouette_plot](https://github.com/user-attachments/assets/4cdf0b49-7925-473e-ac1e-b331663861ff)
![gene_expression_sub_1_umap](https://github.com/user-attachments/assets/bf7d8bfc-966a-4c43-b055-6d230464e6e2)





