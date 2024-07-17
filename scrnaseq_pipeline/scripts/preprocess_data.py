import scanpy as sc
import pandas as pd

def preprocess_data(inputs, output):
    raw_data, metadata_file, cluster_file = inputs.split()
    adata = sc.read_10x_mtx(raw_data).T
    metadata = pd.read_csv(metadata_file, sep='\t', header=None)
    clusters = pd.read_csv(cluster_file, sep='\t', header=None)

    adata.obs['barcode'] = metadata[0].values
    adata.var['gene'] = clusters[1].values
    adata.var_names = adata.var['gene']
    adata.obs_names = adata.obs['barcode']

    adata.write(output)

if __name__ == "__main__":
    import sys
    preprocess_data(" ".join(sys.argv[1:-1]), sys.argv[-1])
