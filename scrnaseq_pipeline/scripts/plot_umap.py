import scanpy as sc

def plot_umap(input, output):
    adata = sc.read_h5ad(input)
    sc.pp.neighbors(adata, n_neighbors=10)
    sc.tl.umap(adata)
    sc.pl.umap(adata, save=output)

if __name__ == "__main__":
    import sys
    plot_umap(sys.argv[1], sys.argv[2])
