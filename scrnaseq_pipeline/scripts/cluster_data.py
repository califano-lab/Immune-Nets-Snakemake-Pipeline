import scanpy as sc

def cluster_data(input, output):
    adata = sc.read_h5ad(input)
    sc.tl.leiden(adata, resolution=0.1)
    adata.write(output)

if __name__ == "__main__":
    import sys
    cluster_data(sys.argv[1], sys.argv[2])
