# Save raw mouse data as Seurat
library(dplyr)
library(Seurat)
library(SeuratDisk)
library(fs) # file system interactions

rm(list = ls())
gc()
# Where to read from and save to 
#datadir <- "~/Dropbox/scRNAseq_data/Bakken21/"
#savedir <- "./data/Bakken21/"
#counts.filename <- "matrix.csv"
#filename.meta <- "metadata.csv"

# Load the data
#counts <- read.table(paste0(datadir, counts.filename), header = TRUE, sep=",") 
#meta.data <- read.csv(paste0(datadir, filename.meta), header=TRUE, skip=0)
#colnames(meta.data)
#tell r to load as integer
datadir <- "~/Dropbox/tPC_data/"
savedir <- "./data/Bakken21/"

setwd(datadir)
Convert(paste0(datadir, "bakken_raw.h5ad"), dest = "h5seurat", overwrite = TRUE)
# wtf: Warning: Unknown file type: h5ad
# Check
human <- LoadH5Seurat(paste0(datadir, "bakken_raw.h5seurat"))

file_move(paste0(datadir, "bakken_raw.h5ad"), 
          paste0(savedir, "bakken.h5ad"))
setwd(datadir)
Convert("bakken_raw.h5ad", dest = "h5seurat", overwrite = TRUE)
setwd("~/Dropbox/elfn1_evolution/")
human <- LoadH5Seurat(paste0(savedir, "bakken.h5seurat"))
