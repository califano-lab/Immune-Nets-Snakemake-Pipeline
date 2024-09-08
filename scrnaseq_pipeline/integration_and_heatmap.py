import scanpy as sc
import pyviper
import matplotlib.pyplot as plt
import numpy as np

def generate_heatmap(data_path, output_path, N=50):
    # Ensure all previous plot windows are closed
    plt.close('all')
    
    data = sc.read_h5ad(data_path)
    
    # Process data to ensure there are no invalid values
    data.X = np.where(data.X <= 0, np.min(data.X[data.X > 0]), data.X)
    
    try:
        # 选择前 N 个最活跃的蛋白质
        protein_set = data.var_names[:N]
        
        # 创建热图
        pyviper.pl.heatmap(data, var_names=protein_set, groupby="leiden", vcenter=0, cmap="RdBu_r", swap_axes=True, show=False)
        
        # 保存图像
        plt.savefig(output_path, bbox_inches="tight")
        plt.close()
        
    except Exception as e:
        print(f"Failed to generate heatmap due to: {e}")

if __name__ == "__main__":
    generate_heatmap(snakemake.input[0], snakemake.output[0])

