import scanpy as sc

def qc(input, output):
    adata = sc.read_h5ad(input)
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    sc.pl.highest_expr_genes(adata, n_top=20, show=False)
    adata.write(output.replace('.html', '.h5ad'))

if __name__ == "__main__":
    import sys
    qc(sys.argv[1], sys.argv[2])
