# Save raw mouse data as Seurat
library(dplyr)
library(Seurat)
library(SeuratDisk)

savedir <- "./data/Bakken21/"
datadir <- "~/Dropbox/transcriptomic-axis/data/csvs_from_ad/Bakken21/"
metadata <- read.csv(paste0(datadir, "obs.csv"), header = TRUE, sep=',')
genenames <- read.csv(paste0(datadir, "var.csv"), header=FALSE, sep=',', skip = 1)
counts <- read.csv(paste0(datadir, "X.csv"), header=FALSE, sep=',')

# Annotate for Seurat object
counts <- t(counts)
rownames(counts) <- genenames$V1
colnames(counts) <- metadata$sample_name
rownames(metadata) <- metadata$sample_name

# Make it a Seurat object
human <- CreateSeuratObject(counts = counts, project = "human", 
                           meta.data = metadata, 
                           min.cells = 1, min.features = 1)
SaveH5Seurat(human, paste0(savedir, "human"), overwrite=TRUE)
# Also save in scanpy/H5ad format
Convert(paste0(savedir, paste0("human", ".h5seurat")), 
        dest = "h5ad", overwrite=TRUE)
