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
    for dataset in ['tosches', 'colquitt', 'tasic', 'bakken']:
        print(f"Processing {dataset}")
        adata = ad.read_h5ad(f"./data/{dataset}_raw.h5ad")
        # Focus on 5 canonical mammalian classes. 
        keep_subclasses = ['Pvalb', 'Sst', 'Lamp5', 'Sncg', 'Vip']
        adata.obs['subclass'] = adata.obs['subclass'].astype("category")
        # Merge Serpinf1 with Vip
        if 'Serpinf1' in adata.obs['subclass'].cat.categories:
            adata.obs['subclass'] = adata.obs['subclass'].cat.remove_categories(['Serpinf1'])
            adata.obs.loc[adata.obs['subclass'].isna(),'subclass'] = 'Vip'
        # Remove Meis2
        extra_categories = [subclass for subclass in adata.obs['subclass'].cat.categories if subclass not in keep_subclasses]
        #for subclass in extra_categories:
        adata = adata[~adata.obs['subclass'].isin(extra_categories)]
        missing_categories = [subclass for subclass in keep_subclasses if subclass not in adata.obs['subclass'].cat.categories]
        adata.obs['subclass'] = adata.obs['subclass'].cat.add_categories(missing_categories)
        print(adata.obs['subclass'].cat.categories)
        #if dataset in ['tosches', 'colquitt']:
        # Order by Elfn1 expression in mice
        adata = adata[adata.obs['subclass']!='Sncg']
        adata.obs['subclass'] = adata.obs['subclass'].cat.reorder_categories(['Sst', 'Vip', 'Pvalb', 'Lamp5'])
        #else:
        #    adata.obs['subclass'] = adata.obs['subclass'].cat.reorder_categories(['Sst', 'Vip', 'Sncg', 'Pvalb', 'Lamp5'])

        # Normalize
        sc.pp.normalize_total(adata, 1e4)
        sc.pp.log1p(adata)
        # Violin plot Elfn1
        fig, ax = plt.subplots(figsize=(9,4))
        sns.despine()
        sc.pl.violin(adata, "Elfn1", "subclass", ax=ax, ylabel = "Elfn1 (log CP10K)", xlabel = "Subclass",
                    save=f"_{dataset}.png", show=False)

        # PCA on HVGs
        sc.pp.highly_variable_genes(adata)
        sc.pp.pca(adata, random_state = seed)
        # UMAP
        sc.pp.neighbors(adata, random_state = seed)
        sc.tl.umap(adata, random_state = seed)

        # Visualize UMAP
        sc.pl.umap(adata, color='subclass', legend_loc='on data', show=False,
                frameon=False, title='Subclass', save=f"_subclass_{dataset}.png")
        sc.pl.umap(adata, color='Elfn1', legend_loc='on data', show = False,
                frameon=False, title='Elfn1', save=f"_elfn1_{dataset}.png")


if __name__ == "__main__":
    main()

