# Save raw bird data as Seurat
library(dplyr)
library(Seurat)
library(SeuratDisk)
rm(list=ls())

# Where to load from and save to
projectdir <- "/home/joram/Dropbox/elfn1_evolution"
if (getwd() != projectdir){
  setwd(projectdir)  
}
source("./scripts/preprocessing_utils.R") # function recode_genes
datadir <- "./data/raw/colquitt/"
savedir <- datadir

# Load the data: single csv with meta data and counts
counts_and_metadata <- read.table(paste0(datadir, "HVC_RA_RNA_counts.csv"), header=TRUE, sep=",")
counts <- counts_and_metadata %>% select(!c(cell,UMAP_1, UMAP_2, X, nCount_RNA,
                                            nFeature_RNA, percent.mito,
                             species, position2, cluster_int_sub2, orig.ident))
meta.data <- counts_and_metadata %>% select(c(cell,UMAP_1, UMAP_2, nCount_RNA, nFeature_RNA, percent.mito,
                             species, position2, cluster_int_sub2, orig.ident))

# Prep data for seurat object. 
counts <- t(counts) # cells x genes -> genes x cells
colnames(counts) <- meta.data$cell # annotate with dim names
rownames(meta.data) <- meta.data$cell
rownames(counts) <- recode_genes(rownames(counts))
meta.data <-meta.data %>% select(!cell) # remove cell ident - not necessary
# Create Seurat
bird <- CreateSeuratObject(counts = counts, project = "bird", 
                           meta.data = meta.data, 
                           min.cells = 1, min.features = 1)
bird # 18342 features across 29497 samples within 1 assay 

SaveH5Seurat(bird, paste0(savedir, "bird"), overwrite=TRUE)
# Also save in scanpy/H5ad format
Convert(paste0(savedir, paste0("bird", ".h5seurat")), 
        dest = "h5ad", overwrite=TRUE)
