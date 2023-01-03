# Visualize integrated datasets
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
import scanpy as sc
from matplotlib.colors import ListedColormap

def main():
    sns.set_context("poster")
    sns.set_palette("Set2")
    data_dir = "./data/anndata/"
    fig_dir = "./figures/"

    for celltype in ['GABAergic', 'Glutamatergic']:
        adata = ad.read_h5ad(f"{data_dir}mouse_bird_{celltype}_integrated_SCT.h5ad")
        fig, ax = plt.subplots()
        sc.pl.umap(adata, color='organism', frameon=False, title=f'{celltype} neurons', 
            ax=ax, legend_loc='None', show=False)
        fig.savefig(f"{fig_dir}integrated_umap_{celltype}_SCT", dpi=300)


if __name__ == "__main__":
    main()

