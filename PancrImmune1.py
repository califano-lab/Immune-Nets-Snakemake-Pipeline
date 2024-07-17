#!/usr/bin/env python
# coding: utf-8

# In[14]:


get_ipython().system('pip install viper-in-python scanpy anndata pandas numpy')


# In[13]:


# 设置文件路径
data_location = "/Users/lzy/Desktop/"
matrix_file = data_location + "expression_matrix.mtx.gz"
metadata_file = data_location + "PancrImmune1-metadata.tsv"
cluster_file = data_location + "PancrImmune1-cluster.tsv"

# 读取并处理 metadata 文件
metadata = pd.read_csv(metadata_file, sep='\t')

# 读取并处理 cluster 文件
cluster = pd.read_csv(cluster_file, header=None, sep='\t')

# 解压缩并读取 matrix 文件
temp_matrix_file = "temp_matrix.mtx"
with gzip.open(matrix_file, 'rt') as f:
    lines = f.readlines()
    with open(temp_matrix_file, 'wt') as temp_f:
        temp_f.writelines(lines[1:])  # 跳过第一行标题

# 读取基因表达签名
gene_expr_signature = sc.read_mtx(temp_matrix_file).T  # 转置矩阵使细胞为行，基因为列

# 获取 matrix 文件中的条形码（假设条形码是行名称）
matrix_barcodes = metadata.iloc[:, 0].values

# 检查 metadata 和 matrix 的条形码是否匹配
metadata_barcodes = metadata.iloc[:, 0].values
if len(metadata_barcodes) != len(set(metadata_barcodes)):
    print("Metadata 中存在重复条形码")
else:
    print("Metadata 中没有重复条形码")

# 找出 matrix 文件中不存在的多余条形码
extra_barcodes = set(metadata_barcodes) - set(matrix_barcodes)

# 删除多余条形码的行
if extra_barcodes:
    print(f"找到多余的条形码: {extra_barcodes}")
    metadata = metadata[~metadata.iloc[:, 0].isin(extra_barcodes)]
    print(f"删除多余条形码后的 metadata 行数: {len(metadata)}")
else:
    print("没有找到多余的条形码")

# 确保 metadata 和表达矩阵行数一致
print(f"Metadata 行数: {len(metadata)}")
print(f"Cluster 文件行数: {len(cluster)}")
print(f"表达矩阵行数: {gene_expr_signature.shape[0]}")
print(f"表达矩阵列数: {gene_expr_signature.shape[1]}")

assert len(metadata) == gene_expr_signature.shape[0], "行数不一致，请检查数据"

# 添加条形码和基因信息
gene_expr_signature.obs['barcode'] = metadata.iloc[:, 0].values  # 假设条形码在第一列
gene_expr_signature.var['gene'] = cluster[1].values
gene_expr_signature.var_names = gene_expr_signature.var['gene']

# 保存为 h5ad 文件
output_file = data_location + "gene_expression_signature.h5ad"
gene_expr_signature.write(output_file)

# 删除临时文件
os.remove(temp_matrix_file)

print(f"处理完成，文件已保存为: {output_file}")


# In[18]:


import pandas as pd

# 设置数据路径
data_location = "/Users/lzy/Desktop/"
metadata_file = data_location + "PancrImmune1-metadata.tsv"

# 读取 metadata 文件
metadata = pd.read_csv(metadata_file, sep='\t', header=None)

# 打印原始 metadata 行数
print(f"Original metadata rows: {len(metadata)}")

# 删除多余的一行（假设多余的行在开头或结尾）
metadata = metadata.drop(metadata.index[0])  # 删除第一行，如果是最后一行，则使用 metadata.drop(metadata.index[-1])

# 打印删除后的 metadata 行数
print(f"Modified metadata rows: {len(metadata)}")

# 保存修改后的 metadata 文件
metadata.to_csv("PancrImmune1-metadata_modified.tsv", sep='\t', index=False, header=False)


# In[1]:


import pandas as pd
import os

# Set data paths
data_location = "/Users/lzy/Desktop/"
metadata_file = os.path.join(data_location, "PancrImmune1-metadata.tsv")
modified_metadata_file = os.path.join(data_location, "PancrImmune1-metadata_modified.tsv")

# Check if the file exists
def check_file(file_path):
    if os.path.exists(file_path):
        print(f"{file_path} exists.")
    else:
        print(f"Error: {file_path} does not exist.")

# Check the original file
check_file(metadata_file)

# Read the metadata file
metadata = pd.read_csv(metadata_file, sep='\t', header=None, low_memory=False)
print(f"Original metadata rows: {len(metadata)}")

# Remove the extra row (assuming the extra row is at the beginning or end)
metadata = metadata.drop(metadata.index[0])  # Remove the first row; if it's the last row, use metadata.drop(metadata.index[-1])
print(f"Modified metadata rows: {len(metadata)}")

# Save the modified metadata file
metadata.to_csv(modified_metadata_file, sep='\t', index=False, header=False)

# Check the modified file
check_file(modified_metadata_file)


# In[4]:


import scanpy as sc
import pandas as pd

# 设置数据路径
data_location = "/Users/lzy/Desktop/"
matrix_file = data_location + "expression_matrix.mtx.gz"
metadata_file = data_location + "PancrImmune1-metadata_modified.tsv"
cluster_file = data_location + "PancrImmune1-cluster.tsv"

# 读取基因表达签名
gene_expr_signature = sc.read_mtx(matrix_file).T  # 转置矩阵使细胞为行，基因为列

# 尝试逐块读取数据文件
try:
    barcodes = pd.read_csv(metadata_file, header=None, sep='\t')
    features = pd.read_csv(cluster_file, header=None, sep='\t')
    
    # 打印维度信息进行检查
    print(f"Dimensions of gene_expr_signature.var: {gene_expr_signature.var.shape}")
    print(f"Dimensions of features: {features.shape}")
    
    # 确保 features 的行数与 gene_expr_signature.var 的行数匹配
    if len(features) != gene_expr_signature.shape[1]:
        raise ValueError("The number of genes in the features file does not match the number of genes in the expression matrix.")
    
    # 添加条形码和基因信息
    gene_expr_signature.obs['barcode'] = barcodes[0].values
    gene_expr_signature.var['gene'] = features[1].values
    gene_expr_signature.var_names = gene_expr_signature.var['gene']
    gene_expr_signature.obs_names = gene_expr_signature.obs['barcode']

    print(gene_expr_signature)

    # 计算邻近细胞
    sc.pp.neighbors(gene_expr_signature, n_neighbors=10)

    # 计算UMAP
    sc.tl.umap(gene_expr_signature)

    # 绘制UMAP图
    sc.pl.umap(gene_expr_signature, color=['donor', 'sex', 'cell_type__ontology_label'])

    # 使用Leiden算法聚类
    sc.tl.leiden(gene_expr_signature, resolution=0.1)

    # 绘制聚类结果
    sc.pl.umap(gene_expr_signature, color='leiden')
    
except FileNotFoundError as e:
    print(f"Error: {e}")
    print("请检查文件路径和文件名是否正确，并确保文件存在。")
except ValueError as e:
    print(f"ValueError: {e}")


# In[2]:


import pyviper
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import random
import warnings

warnings.filterwarnings("ignore", message=".*The 'nopython' keyword.*")  # for jit decorator issue with sc.pp.neighbors

# 设置数据路径
data_location = "/Users/lzy/Desktop/"
matrix_file = data_location + "expression_matrix.mtx.gz"
metadata_file = data_location + "PancrImmune1-metadata_modified.tsv"  # 使用修改后的metadata文件
cluster_file = data_location + "PancrImmune1-cluster.tsv"


# 读取基因表达签名
gene_expr_signature = sc.read_mtx(matrix_file).T  # 转置矩阵使细胞为行，基因为列
barcodes = pd.read_csv(metadata_file, header=None, sep='\t')
features = pd.read_csv(cluster_file, header=None, sep='\t')

# 添加条形码和基因信息
gene_expr_signature.obs['barcode'] = barcodes[0].values
gene_expr_signature.var['gene'] = features[1].values
gene_expr_signature.var_names = gene_expr_signature.var['gene']
gene_expr_signature.obs_names = gene_expr_signature.obs['barcode']

print(gene_expr_signature)

# 计算邻近细胞
sc.pp.neighbors(gene_expr_signature, n_neighbors=10)

# 计算UMAP
sc.tl.umap(gene_expr_signature)

# 绘制UMAP图
sc.pl.umap(gene_expr_signature, color=['donor', 'sex', 'cell_type__ontology_label'])

# 使用Leiden算法聚类
sc.tl.leiden(gene_expr_signature, resolution=0.1)

# 绘制聚类结果
sc.pl.umap(gene_expr_signature, color='leiden')


# In[20]:


get_ipython().system('pip install viper-in-python')


# In[21]:


import pyviper
import scanpy as sc
import anndata 
import pandas as pd
import numpy as np
import random
import warnings

warnings.filterwarnings("ignore", message=".*The 'nopython' keyword.*")  # For jit decorator issue with sc.pp.neighbors

