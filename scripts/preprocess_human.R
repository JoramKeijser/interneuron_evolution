# Save human data as Seurat
library(dplyr)
library(Seurat)
library(SeuratDisk)
rm(list=ls())

# Where to load from and save to
projectdir <- "/home/joram/Dropbox/elfn1_evolution"
setwd(projectdir)  
source("./src/preprocessing_utils.R") # function recode_genes & save dirs
datadir <- "./data/raw/bakken/"

nrows = -1 
metadata <- read.csv(paste0(datadir, "metadata.csv"), header = TRUE, sep=',', 
                     nrows = nrows)
counts <- read.csv(paste0(datadir, "matrix.csv"), header=TRUE, sep=',', nrows=nrows)
# sample name as rowname
rownames(counts) <- counts$sample_name
counts <- counts %>% select(!sample_name)
# Annotate for Seurat object
counts <- t(counts)
colnames(counts) <- metadata$sample_name
rownames(metadata) <- metadata$sample_name
# Gene names in "Elfn1" format (mouse convention)
rownames(counts) <- recode_genes(rownames(counts))


# Make it a Seurat object
human <- CreateSeuratObject(counts = counts, project = "human", 
                           meta.data = metadata, 
                           min.cells = 1, min.features = 1)

# Save it, also has h5ad for python use
SaveH5Seurat(mouse, paste0(savedir_seurat, "human"), overwrite=TRUE)
# Also save in scanpy/H5ad format
Convert(paste0(savedir_seurat, "human.h5seurat"), 
        dest = paste0(savedir_anndata, "human.h5ad"), overwrite=TRUE)

