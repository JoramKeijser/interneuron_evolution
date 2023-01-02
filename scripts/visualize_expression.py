import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
import scanpy as sc
from matplotlib.colors import ListedColormap
sns.set_context("poster")
sns.set_palette("colorblind")
sc.settings.figdir = "./figures/"

def main():
    seed = 1029
    data_dir = "../tPC_data/"
    # Upper lim. mapped to color. Higher expressions are clipped. 
    max_expr = {'tosches': 2.2, 'colquitt': 4.2, "tasic": 1.75, "bakken": 1.25}
    for dataset in ['tasic', 'bakken', 'tosches', 'colquitt']:
        print(f"Processing {dataset}")
        adata = ad.read_h5ad(f"{data_dir}{dataset}_raw.h5ad")
        # Focus on 5 canonical mammalian classes. 
        keep_subclasses = ['Pvalb', 'Sst', 'Lamp5', 'Sncg', 'Vip']
        adata.obs['subclass'] = adata.obs['subclass'].astype("category")
        # Merge Serpinf1 with Vip
        if 'Serpinf1' in adata.obs['subclass'].cat.categories:
            adata.obs['subclass'] = adata.obs['subclass'].cat.remove_categories(['Serpinf1'])
            adata.obs.loc[adata.obs['subclass'].isna(),'subclass'] = 'Vip'
        # Only keep canonical subclasses (remove e.g. Meis2)
        extra_categories = [subclass for subclass in adata.obs['subclass'].cat.categories if subclass not in keep_subclasses]
        adata = adata[~adata.obs['subclass'].isin(extra_categories)]
        missing_categories = [subclass for subclass in keep_subclasses if subclass not in adata.obs['subclass'].cat.categories]
        adata.obs['subclass'] = adata.obs['subclass'].cat.add_categories(missing_categories)
        adata = adata[adata.obs['subclass']!='Sncg']
        # Order by Elfn1 expression in mice
        adata.obs['subclass'] = adata.obs['subclass'].cat.reorder_categories(['Sst', 'Vip', 'Pvalb', 'Lamp5'])
        #else:
        #    adata.obs['subclass'] = adata.obs['subclass'].cat.reorder_categories(['Sst', 'Vip', 'Sncg', 'Pvalb', 'Lamp5'])

        # Normalize to log CP10K
        sc.pp.normalize_total(adata, 1e4)
        sc.pp.log1p(adata)
        # Violin plots Elfn1 amd Cbln4
        for gene in ['Elfn1', 'Cbln4']:
            fig, ax = plt.subplots(figsize=(9,4))
            sns.despine()
            sc.pl.violin(adata, gene, "subclass", ax=ax, ylabel = f"{gene} (log CP10K)", xlabel = "Subclass",
                        save=f"_{gene}_{dataset}.png", show=False)

        # PCA on HVGs
        sc.pp.highly_variable_genes(adata)
        sc.pp.pca(adata, random_state = seed)
        # UMAP
        sc.pp.neighbors(adata, random_state = seed)
        sc.tl.umap(adata, random_state = seed)

        # Visualize UMAP
        sc.pl.umap(adata, color='subclass', legend_loc='on data', show=False,
                frameon=False, title='Subclass', save=f"_subclass_{dataset}.png")
        for gene in ['Elfn1', 'Cbln4', 'Calb2']:
            sc.pl.umap(adata, color=gene, legend_loc='on data', show = False,
                frameon=False, title=gene, save=f"_{gene}_{dataset}.png", vmax=max_expr[dataset])


if __name__ == "__main__":
    main()

