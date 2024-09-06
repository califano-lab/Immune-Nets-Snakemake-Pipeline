# Project Title: scRNA-seq Data Analysis Pipeline with PyVIPER

This repository contains a pipeline for analyzing single-cell RNA-seq (scRNA-seq) data. The analysis includes preprocessing, quality control, normalization, gene regulatory network inference using ARACNe3, and regulon activity analysis using PyViper.

## Overview
The goal of this pipeline is to process gene expression datasets and analyze them to understand regulatory networks and protein activity in single cells. This process involves multiple steps, including data preparation, ARACNe3 network generation, VIPER analysis, PCA, clustering, and visualization.

## Introduction of Snakemake

Snakemake is a workflow management system that allows you to create reproducible and scalable data analyses. Inspired by the Python programming language and the "Make" utility commonly used for build automation, Snakemake enables you to define complex workflows using simple rules that specify how input files are converted into output files.

## Pipeline Overview

```
# Read input and output folder from config or command-line arguments
input_folder = config.get("input_folder", "data")  # Default is "data" if not provided
output_folder = config.get("output_folder", "results")  # Default is "results" if not provided
​
# Detect all gene expression files in the input folder dynamically
import os
import glob
​
gene_expr_files = sorted(glob.glob(f"{input_folder}/*.tsv"))
datasets = [os.path.splitext(os.path.basename(f))[0] for f in gene_expr_files]
​
rule all:
    input:
        expand(f"{output_folder}/{{dataset}}_heatmap.png", dataset=datasets),
        expand(f"{output_folder}/{{dataset}}_umap.png", dataset=datasets)
​
# Step 1: Convert and rename files if necessary
rule rename_files:
    input:
        gene_expr_dir=input_folder
    output:
        "renamed.complete"
    script:
        "scripts/rename_files.py"
​
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
​
# Step 3: Convert ARACNe3 output to formatted networks
rule format_network:
    input:
        aracne_output=f"{output_folder}/{{dataset}}_consolidated-net_defaultid/consolidated-net_defaultid.tsv"
    output:
        formatted_network=f"{output_folder}/{{dataset}}_consolidated-net_defaultid_formatted_network.tsv"
    script:
        "scripts/format_network.py"
​
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
​
# Step 5: VIPER analysis and PCA
rule viper_and_pca_analysis:
    input:
        processed_expr=f"{output_folder}/{{dataset}}_processed_expr.h5ad",
        processed_net=f"{output_folder}/{{dataset}}_processed_net.pkl"
    output:
        prot_act_pca=f"{output_folder}/{{dataset}}_prot_act_pca.h5ad"
    script:
        "scripts/viper_and_pca_analysis.py"
​
# Step 6: Clustering and UMAP analysis
rule clustering_and_umap:
    input:
        prot_act_pca=f"{output_folder}/{{dataset}}_prot_act_pca.h5ad"
    output:
        umap_data=f"{output_folder}/{{dataset}}_umap_data.h5ad",
        umap_plot=f"{output_folder}/{{dataset}}_umap.png"
    script:
        "scripts/clustering_and_umap.py"
​
# Step 7: Generate Heatmap
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

Outputs:
· Heatmap images for each dataset (```results/{dataset}_heatmap.png```).
· UMAP plot images for each dataset (```results/{dataset}_umap.png```).

Functionality: Acts as the endpoint to verify that all datasets have been processed through the pipeline.

2.rename_files:

Purpose: This optional step standardizes the names and formats of gene expression files in the data directory to ensure consistency.

Inputs:
· The directory containing all the gene expression files (```gene_expr_dir="data"```).

Outputs:

A completion flag file (```renamed.complete```) to indicate that all files have been renamed and formatted correctly.

Script: ```scripts/rename_files.py``` is used to perform the renaming and formatting operations.

3.run_aracne3:

Purpose: Executes ARACNe3, a tool for inferring gene regulatory networks, using gene expression data and a combined regulators file.

Inputs:

· Gene expression file (```data/{dataset}.tsv```).
· Combined regulators file (```combined_regulators.txt```).

Outputs:

A directory containing the ARACNe3 consolidated network output (```results/{dataset}_consolidated-net_defaultid```).

Shell Command: Runs the ARACNe3 executable with specified parameters to generate the network.

4.format_network:

Purpose: Converts the output from ARACNe3 into a format that is compatible with PyViper for subsequent analysis.

Inputs:
· ARACNe3 output file (```results/{dataset}_consolidated-net_defaultid/consolidated-net_defaultid.tsv```).

Outputs:
A formatted network file (```results/{dataset}_consolidated-net_defaultid_formatted_network.tsv```).

Script: Uses ```scripts/format_network.py``` to format the network data into a usable structure for PyViper.

5.load_and_preprocess:

Purpose: Loads gene expression data and formatted network data, and performs preprocessing necessary for downstream VIPER analysis.

Inputs:

Gene expression file (```data/{dataset}.tsv```).
Formatted network file (```results/{dataset}_consolidated-net_defaultid_formatted_network.tsv```).

Outputs:
· Preprocessed gene expression data (```results/{dataset}_processed_expr.h5ad```).
· Processed network file (```results/{dataset}_processed_net.pkl```).

Script: Executes ```scripts/load_and_preprocess.py``` to preprocess the data.

6. viper_and_pca_analysis:

Purpose: Conducts VIPER analysis to infer protein activity from gene expression data and applies PCA for dimensionality reduction.

Inputs:
Preprocessed gene expression data (```results/{dataset}_processed_expr.h5ad```).
Processed network file (```results/{dataset}_processed_net.pkl```).

Outputs:
Protein activity PCA data (```results/{dataset}_prot_act_pca.h5ad```).
Script: Runs ```scripts/viper_and_pca_analysis.py``` to perform VIPER and PCA analysis..

7. clustering_and_umap:

Purpose: Performs clustering analysis and generates UMAP visualizations to visualize the data in two-dimensional space.

Inputs:

Protein activity PCA data (```results/{dataset}_prot_act_pca.h5ad```)

Outputs:

UMAP data (```results/{dataset}_umap_data.h5ad```).
UMAP plot image (```results/{dataset}_umap.png```).

Script: Executes ```clustering_and_umap.py``` to perform clustering and generate UMAP plots.

8. integration_and_heatmap:

Purpose: Integrates the results from UMAP and clustering steps and generates a heatmap to visualize protein activity across clusters.

Inputs:
UMAP data (```results/{dataset}_umap_data.h5ad```).

Outputs:

Final heatmap image (```results/{dataset}_heatmap.png```).
Script: Runs ```scripts/integration_and_heatmap.py``` to create a heatmap visualization.


## Code Part

## 1. rename files

```
import os
import glob
import pandas as pd

# Define the data directory
data_dir = "/Users/lzy/Desktop/Final/data"

# Get all files in the directory (excluding those already in TSV format)
files = sorted(glob.glob(os.path.join(data_dir, "*")))
tsv_files = [f for f in files if f.endswith('.tsv')]

# Iterate over each file and convert to TSV format if necessary
for file_path in files:
    if file_path in tsv_files:
        print(f"File {file_path} is already in TSV format. Skipping conversion.")
        continue

    # Read the file
    try:
        # Assume the file is in CSV format or another format separated by commas
        df = pd.read_csv(file_path, sep=None, engine='python')
        
        # Get the base name of the file (without the extension)
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        
        # Generate the new TSV file path
        new_name = os.path.join(data_dir, f"{base_name}.tsv")
        
        # Save the data as a TSV file
        df.to_csv(new_name, sep='\t', index=False)
        print(f"Converted {file_path} to {new_name}.")
        
        # Delete the original file
        os.remove(file_path)
        print(f"Deleted original file {file_path}.")
    except Exception as e:
        print(f"Error processing {file_path}: {e}")

print("All files have been successfully converted to TSV format and original files have been deleted.")
```

## 2. format_network
Shell code:
```
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

## 3. load_and_preprocess


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

##  4. Viper and pca analysis

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

##  5. UMAP and Clustering

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

    # Ensure the output directory exists
    os.makedirs(os.path.dirname(output_figure), exist_ok=True)

    # Save the UMAP plot
    sc.pl.umap(data, color='leiden', show=False, save=False)  # Generate the plot without directly saving it
    plt.savefig(output_figure)  # Save using matplotlib to ensure the path is correct
    plt.close()

    # Save the processed data
    data.write_h5ad(output_file)

if __name__ == "__main__":
    perform_clustering_and_umap(
        snakemake.input.prot_act_pca, 
        snakemake.output.umap_data, 
        snakemake.output.umap_plot
    )


```
##  6. Heatmap Generation

```
import scanpy as sc
import pyviper
import matplotlib.pyplot as plt
import numpy as np

def generate_heatmap(data_path, output_path):
    # Ensure all previous plot windows are closed
    plt.close('all')
    
    data = sc.read_h5ad(data_path)
    
    # Process data to ensure there are no invalid values
    data.X = np.where(data.X <= 0, np.min(data.X[data.X > 0]), data.X)

    # Generate and save the heatmap
    try:
        sc.pp.neighbors(data, n_neighbors=10, n_pcs=40)
        sc.tl.leiden(data)
        sc.tl.umap(data)
        sc.tl.rank_genes_groups(data, 'leiden', method='t-test', n_genes=20)
        
        # Generate the heatmap using Scanpy's heatmap plotting function
        sc.pl.heatmap(data, var_names=data.uns['rank_genes_groups']['names']['0'], groupby='leiden', cmap='viridis', show=False)
        
        # Save the image
        plt.savefig(output_path)
        plt.close()  # Close the image to ensure no additional blank plots are generated
        
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

![gene_expression_sub_10_umap](https://github.com/user-attachments/assets/a3b62ed8-57a8-4b0a-b81f-4cb02f9e291e)
![gene_expression_sub_10_heatmap](https://github.com/user-attachments/assets/beec4387-0ce3-4e98-87c7-d240ed13ccd9)
![gene_expression_sub_9_umap](https://github.com/user-attachments/assets/05740084-328b-4fd3-a491-3cedc4947818)
![gene_expression_sub_9_heatmap](https://github.com/user-attachments/assets/b2e89743-4606-404c-960d-5eebcef986d9)
![gene_expression_sub_8_umap](https://github.com/user-attachments/assets/175a4dd1-5a0c-403c-81a6-f0e8f08e6b7e)
![gene_expression_sub_8_heatmap](https://github.com/user-attachments/assets/8689d8d6-e006-45cf-9249-f0c771b0eace)
![gene_expression_sub_7_umap](https://github.com/user-attachments/assets/f364b470-6511-45c8-b87e-ac123c86e794)
![gene_expression_sub_7_heatmap](https://github.com/user-attachments/assets/e65449a8-3191-4e44-b0da-05158304b8cf)
![gene_expression_sub_6_umap](https://github.com/user-attachments/assets/bc01af72-158e-4f79-82cf-994114d9a2ef)
![gene_expression_sub_6_heatmap](https://github.com/user-attachments/assets/9d3c6d4b-b56d-4791-870d-6402fe6e6bfe)
![gene_expression_sub_5_umap](https://github.com/user-attachments/assets/3b2a9286-3ac8-437d-b179-25732e8b6a1b)
![gene_expression_sub_5_heatmap](https://github.com/user-attachments/assets/c244688e-a18e-429c-bd7f-30615060d655)
![gene_expression_sub_4_umap](https://github.com/user-attachments/assets/57904012-f296-4a0a-acc5-d4af54d0d96e)
![gene_expression_sub_4_heatmap](https://github.com/user-attachments/assets/f1a0f520-51a2-4b9c-aaf2-dc3c71ecdb1a)
![gene_expression_sub_3_umap](https://github.com/user-attachments/assets/2c4d882a-3fd0-4004-97c6-e2058ca6fe2d)
![gene_expression_sub_3_heatmap](https://github.com/user-attachments/assets/bf523060-c414-4169-8e0a-5b8bcdaf5ae7)
![gene_expression_sub_2_umap](https://github.com/user-attachments/assets/a146ea29-7f3c-4460-8824-75bdb5772742)
![gene_expression_sub_2_heatmap](https://github.com/user-attachments/assets/2f54a6f3-95d3-4c3e-bb7f-e5c920923fef)
![gene_expression_sub_1_umap](https://github.com/user-attachments/assets/45c4497b-414a-457d-bd7b-34b8c0e13cb9)
![gene_expression_sub_1_heatmap](https://github.com/user-attachments/assets/a802e9dd-13c3-42d3-956e-cab610f87b1c)







