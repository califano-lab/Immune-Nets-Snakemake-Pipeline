import pyviper
import anndata as ad
import pandas as pd

def run_viper(input, output):
    aracne_results = pd.read_csv(input, sep='\t')
    adata = ad.AnnData()
    adata.write(output)

if __name__ == "__main__":
    import sys
    run_viper(sys.argv[1], sys.argv[2])

