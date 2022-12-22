# Save raw bird data as Seurat
library(dplyr)
library(Seurat)
library(SeuratDisk)

datadir <- "~/Dropbox/scRNAseq_data/Colquitt21/"
savedir <- "./data/Colquitt21/"
bird <- read.table(paste0(datadir, "HVC_RA_RNA_counts.csv"), header=TRUE, sep=",")
counts <- bird %>% select(!c(cell,UMAP_1, UMAP_2, X, nCount_RNA, nFeature_RNA, percent.mito,
                             species, position2, cluster_int_sub2, orig.ident))
meta.data <- bird %>% select(c(cell,UMAP_1, UMAP_2, nCount_RNA, nFeature_RNA, percent.mito,
                             species, position2, cluster_int_sub2, orig.ident))

# Prep data for seurat object. 
# "Row names in the metadata need to match the column names of the counts matrix"
# Metadata cells x vars. Counts: genes x cells
counts <- t(counts) # cells x genes -> genes x cells
colnames(counts) <- meta.data$cell
rownames(meta.data) <- meta.data$cell
meta.data <-meta.data %>% select(!cell)
# Create Seurat
bird <- CreateSeuratObject(counts = counts, project = "bird", 
                           meta.data = meta.data, 
                           min.cells = 1, min.features = 1)
bird # 18342 features across 29497 samples within 1 assay 
SaveH5Seurat(bird, paste0(savedir, "bird"), overwrite=TRUE)
# Also save in scanpy/H5ad format
Convert(paste0(savedir, paste0("bird", ".h5seurat")), 
        dest = "h5ad", overwrite=TRUE)
