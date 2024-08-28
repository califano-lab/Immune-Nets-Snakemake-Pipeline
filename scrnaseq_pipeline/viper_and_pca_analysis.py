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
