# Save raw bird data as Seurat
library(dplyr)
library(Seurat)
library(SeuratDisk)
source("./src/preprocessing_utils.R") # function recode_genes & save dirs
datadir <- "./data/raw/tasic/"

# The files to load
filename <- "mouse_VISp_2018-06-14_exon-matrix.csv"
filename.meta <- "mouse_VISp_2018-06-14_samples-columns.csv"
filename.genes <- "mouse_VISp_2018-06-14_genes-rows.csv"
# Load the data
counts <- read.table(paste0(datadir, filename), header = TRUE, sep=",") 
gene.info <- read.table(paste0(datadir, filename.genes), header=TRUE, sep=',')
meta.data <- read.table(paste0(datadir, filename.meta), header=TRUE, sep=',')

# samples as row names. gene x sample
rownames(counts) <- counts$X
counts <- counts %>% select(!X) #genes x samples = 45768 x 15413
rownames(meta.data) <- colnames(counts) # Meta / obs data
rownames(counts) <- gene.info$gene_symbol # Gene names


# Create Seurat
# "Row names in the metadata need to match the column names of the counts matrix"
mouse <- CreateSeuratObject(counts = counts, project = "mouse", 
                            meta.data = meta.data, 
                            min.cells = 1, min.features = 1)
mouse # 42094 features across 15413 samples within 1 assay 

# Save it, also has h5ad for python use
SaveH5Seurat(mouse, paste0(savedir_seurat, "mouse"), overwrite=TRUE)
# Also save in scanpy/H5ad format
Convert(paste0(savedir_seurat, "mouse.h5seurat"), 
        dest = paste0(savedir_anndata, "mouse.h5ad"), overwrite=TRUE)





