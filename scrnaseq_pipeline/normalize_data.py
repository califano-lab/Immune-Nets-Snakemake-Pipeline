import scanpy as sc

data_path = "/Users/lzy/Desktop/"
adata = sc.read(data_path + "qc_data.h5ad")

# Normalization and transformation
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var['highly_variable']]
sc.pp.regress_out(adata, ['total_counts', 'percent_mito'])
sc.pp.scale(adata, max_value=10)

# Save normalized data
adata.write(data_path + "normalized_data.h5ad")