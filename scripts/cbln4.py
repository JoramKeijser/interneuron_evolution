import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
import scanpy as sc
from matplotlib.colors import ListedColormap
sns.set_context("poster")
sns.set_palette("colorblind")
sc.settings.figdir = "./figures/"

seed = 1029 # set randomness for dim. reduction and clustering
#data_dir = "../data/anndata/" 
#dataset='mouse'
#adata = ad.read_h5ad(f"{data_dir}{dataset}.h5ad") #  15413 × 42094
data_dir = "../transcriptomic-axes/results/anndata/"
dataset='tasic'
adata = ad.read_h5ad(f"{data_dir}{dataset}_raw.h5ad")
print(adata)
print(adata.var_names)

print("Preprocess data")
sc.pp.normalize_total(adata, 1e4)
sc.pp.log1p(adata)
# Select Sst cells w/o Chodl 
sst = adata[(adata.obs['subclass'] == "Sst")]
sst = sst[sst[:,"Chodl"].X < 0.5]

print("Dimensionality reduction")
sc.pp.highly_variable_genes(sst)
sc.pp.pca(sst, random_state = seed)
# UMAP
sc.pp.neighbors(sst, random_state = seed)
sc.tl.umap(sst, random_state = seed)
# Clustering
sc.tl.leiden(sst)

# Plot
print("Plotting")
sc.pl.umap(sst, color='leiden', legend_loc='on data',
        frameon=False, title='Clusters', save=f"_sst_cluster.png")
sc.pl.umap(sst, color='brain_subregion', 
        frameon=False, title='Layer', save=f"_sst_layer.png")
gene = 'Cbln4'
sc.pl.umap(sst, color=gene, 
        frameon=False, title=gene, vmax=1.75, save=f"_sst_{gene}.png")

# Marker genes
sc.pl.umap(sst, color=['Tac1', 'Calb2', 'Etv1'],
    frameon=False, vmax=1.75, save=f"_sst_markers.png")