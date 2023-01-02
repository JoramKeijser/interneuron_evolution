# Preprocess turtle data:
# raw to h5ad
import numpy as np 
import anndata as ad
import pandas as pd


def main():
    datadir = "./data/raw/tosches/"
    savedir = "./data/anndata/"
    # Load the data
    anndata = ad.read_h5ad(datadir + "turtle.h5ad")
    print(anndata)
    
    return 0


if __name__ == "__main__":
    main()