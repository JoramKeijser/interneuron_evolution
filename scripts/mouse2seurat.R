# Save raw mouse data as Seurat
library(dplyr)
library(Seurat)
library(SeuratDisk)

rm(list = ls())
gc()
# Where to read from and save to 
datadir <- "~/Dropbox/scRNAseq_data/Tasic18/mouse_VISp_gene_expression_matrices_2018-06-14/"
savedir <- "./data/Tasic18/"
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
# Gene names
# Meta / obs data
rownames(meta.data) <- colnames(counts)
# code as upper case for consistency with mouse dataset
rownames(counts) <- lapply(gene.info$gene_symbol, toupper)


# Create Seurat
# "Row names in the metadata need to match the column names of the counts matrix"
mouse <- CreateSeuratObject(counts = counts, project = "mouse", 
                            meta.data = meta.data, 
                            min.cells = 1, min.features = 1)
mouse # 42094 features across 15413 samples within 1 assay 
SaveH5Seurat(mouse, paste0(savedir, "mouse"), overwrite=TRUE)
# Also save in H5ad format for scanpy
Convert(paste0(savedir, paste0("mouse", ".h5seurat")), 
        dest = "h5ad", overwrite=TRUE)





