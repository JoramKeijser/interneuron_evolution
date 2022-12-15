# Save raw bird data as Seurat
library(dplyr)
library(Seurat)
library(SeuratDisk)
SAVE_RAW = FALSE

bird <- read.table("./data/raw/bird/HVC_RA_RNA_counts.csv", header=TRUE, sep=",")
counts <- bird %>% select(!c(cell,UMAP_1, UMAP_2, X, nCount_RNA, nFeature_RNA, percent.mito,
                             species, position2, cluster_int_sub2, orig.ident))
meta.data <- bird %>% select(c(cell,UMAP_1, UMAP_2, nCount_RNA, nFeature_RNA, percent.mito,
                             species, position2, cluster_int_sub2, orig.ident))
if (SAVE_RAW){
  write.table(meta.data, "./data/raw/bird/metadata.csv", sep=",")
  write.table(counts, "./data/raw/bird/counts.csv", sep=",")  
}

# Create Seurat
# "Row names in the metadata need to match the column names of the counts matrix"
# Metadata cells x vars. Counts: genes x cells
counts <- t(counts) # cells x genes -> genes x cells
colnames(counts) <- meta.data$cell

rownames(meta.data) <- meta.data$cell
meta.data <-meta.data %>% select(!cell)
bird <- CreateSeuratObject(counts = counts, project = "bird", 
                           meta.data = meta.data, 
                           min.cells = 1, min.features = 1)
bird # 18342 features across 29497 samples within 1 assay 
SaveH5Seurat(bird, "./data/raw/bird/bird", overwrite=TRUE)

# SCTransform
bird <- SCTransform(bird, vars.to.regress = "percent.mito", 
                    verbose = TRUE, variable.features.n = 3000)
SaveH5Seurat(bird, "./data/raw/bird/birdSCT", overwrite=TRUE)

