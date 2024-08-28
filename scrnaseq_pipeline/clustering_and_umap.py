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

    # Save the UMAP figure
    sc.pl.umap(data, color='leiden', show=False, save=False)  # Generate the plot without directly saving
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
