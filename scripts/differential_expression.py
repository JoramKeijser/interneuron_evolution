import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
import scanpy as sc


def main():
    for dataset in ['bakken', 'tosches', 'colquitt', 'tasic']:
        print(f"Processing {dataset}")
        adata = ad.read_h5ad(f"./data/{dataset}_raw.h5ad")
        # Normalize
        sc.pp.normalize_per_cell(adata, 1e4)
        sc.pp.log1p(adata)
        # Compare Sst with Pvalb
        print("Determining marker genes")
        sc.tl.rank_genes_groups(adata, "subclass", groups=['Sst', 'Pvalb'], 
            method='t-test', n_genes=1000)

        result = adata.uns['rank_genes_groups']
        groups = result['names'].dtype.names
        df = pd.DataFrame({group + '_' + key[:1]: result[key][group]
            for group in groups for key in ['names', 'logfoldchanges']})
        logfc = float(df[df['Sst_n'] == 'Elfn1']['Sst_l'])
        rank  = np.where(df['Sst_n'] == 'Elfn1')[0][0] + 1
        print(f"{dataset}: Elfn1 is {rank} with {logfc:0.2f} logfc")

    return 0


if __name__ == "__main__":
    main()

