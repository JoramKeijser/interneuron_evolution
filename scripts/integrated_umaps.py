import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
import scanpy as sc
from matplotlib.colors import ListedColormap
sns.set_context("poster")
sns.set_palette("Set2")
sc.settings.figdir = "./figures/"

def main():
    
    data_dir = "./data/raw/"
    fig_dir = "./figures/"

    for celltype in ['GABAergic', 'Glutamatergic']:
        adata = ad.read_h5ad(f"{data_dir}mouse_bird_{celltype}_integrated.h5ad")
        fig, ax = plt.subplots()
        sc.pl.umap(adata, color='organism', frameon=False, title=f'{celltype} neurons', 
            ax=ax, legend_loc='None', show=False)
        fig.savefig(f"{fig_dir}integrated_umap_{celltype}", dpi=300)


if __name__ == "__main__":
    main()

