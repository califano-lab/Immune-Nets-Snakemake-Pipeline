import scanpy as sc
import pyviper
import matplotlib.pyplot as plt
import numpy as np

def generate_heatmap(data_path, output_path):
    # Ensure all previous plot windows are closed
    plt.close('all')
    
    data = sc.read_h5ad(data_path)
    
    # Process the data to ensure there are no invalid values
    data.X = np.where(data.X <= 0, np.min(data.X[data.X > 0]), data.X)

    # Generate and save the heatmap
    try:
        sc.pp.neighbors(data, n_neighbors=10, n_pcs=40)
        sc.tl.leiden(data)
        sc.tl.umap(data)
        sc.tl.rank_genes_groups(data, 'leiden', method='t-test', n_genes=20)
        
        # Use scanpy's heatmap plotting function to generate a heatmap
        sc.pl.heatmap(data, var_names=data.uns['rank_genes_groups']['names']['0'], groupby='leiden', cmap='viridis', show=False)
        
        # Save the image
        plt.savefig(output_path)
        plt.close()  # Close the plot to ensure no extra blank images are generated
        
    except Exception as e:
        print(f"Failed to generate heatmap due to: {e}")

if __name__ == "__main__":
    generate_heatmap(snakemake.input[0], snakemake.output[0])
