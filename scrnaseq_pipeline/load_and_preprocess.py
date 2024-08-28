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
