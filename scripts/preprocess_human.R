# Save human data as Seurat
library(dplyr)
library(Seurat)
library(SeuratDisk)
rm(list=ls())

# Where to load from and save to
projectdir <- "/home/joram/Dropbox/elfn1_evolution"
if (getwd() != projectdir){
  setwd(projectdir)  
}
source("./src/preprocessing_utils.R") # function recode_genes
datadir <- "./data/raw/bakken/"
savedir <- datadir
nrows = 100
metadata <- read.csv(paste0(datadir, "metadata.csv"), header = TRUE, sep=',', 
                     nrows = nrows)
genenames <- read.csv(paste0(datadir, "var.csv"), header=FALSE, sep=',', 
                      skip = 1, nrows = nrows)
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

SaveH5Seurat(human, paste0(savedir, "human"), overwrite=TRUE)
# Also save in scanpy/H5ad format
Convert(paste0(savedir, paste0("human", ".h5seurat")), 
          dest = "h5ad", overwrite=TRUE)
