# Save raw bird data as Seurat
# https://satijalab.org/seurat/articles/integration_introduction.html
library(dplyr)
library(Seurat)
library(SeuratDisk)

rm(list = ls())
filename <- "mouse_VISp_2018-06-14_exon-matrix.csv"
filename.meta <- "mouse_VISp_2018-06-14_samples-columns.csv"
filename.genes <- "mouse_VISp_2018-06-14_genes-rows.csv"
datadir <- "./data/raw/mouse/"
counts <- read.table(paste0(datadir, filename), header = TRUE, sep=",") 

# samples as row names. gene x sample
rownames(counts) <- counts$X
counts <- counts %>% select(!X) #genes x samples = 45768 x 15413
# Gene names
gene.info <- read.table(paste0(datadir, filename.genes), header=TRUE, sep=',')
# Meta / obs data
meta.data <- read.table(paste0(datadir, filename.meta), header=TRUE, sep=',')
rownames(meta.data) <- colnames(counts)
# code as upper case for consistency with mouse dataset
rownames(counts) <- lapply(gene.info$gene_symbol, toupper)

# Create Seurat
# "Row names in the metadata need to match the column names of the counts matrix"
mouse <- CreateSeuratObject(counts = counts, project = "mouse", 
                            meta.data = meta.data, 
                            min.cells = 1, min.features = 1)
mouse # 42094 features across 15413 samples within 1 assay 
SaveH5Seurat(mouse, "./data/raw/mouse/mouse", overwrite=TRUE)
# Transform
mouse <- SCTransform(mouse, vars.to.regress = "percent_mt_exon_reads", 
                     verbose=TRUE, variable.features.n = 3000)
SaveH5Seurat(mouse, "./data/raw/mouse/mouseSCT", overwrite=TRUE)

# Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')



