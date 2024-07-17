import anndata as ad
import pandas as pd

def run_aracne(input, output):
    adata = ad.read_h5ad(input)
    pd.DataFrame().to_csv(output, sep='\t')

if __name__ == "__main__":
    import sys
    run_aracne(sys.argv[1], sys.argv[2])
