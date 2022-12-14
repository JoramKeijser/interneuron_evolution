# Integrate datasets
# https://satijalab.org/seurat/articles/integration_introduction.html
library(dplyr)
library(Seurat)
library(SeuratDisk)

bird <- read.table("./data/raw/bird/HVC_RA_RNA_counts.csv", header=TRUE, sep=",")
head(bird)
colnames(bird)
rownames(bird)

counts <- bird %>% select(!c(cell,UMAP_1, UMAP_2, X, nCount_RNA, nFeature_RNA, percent.mito,
                             species, position2, cluster_int_sub2, orig.ident))
meta.data <- bird %>% select(c(cell,UMAP_1, UMAP_2, nCount_RNA, nFeature_RNA, percent.mito,
                             species, position2, cluster_int_sub2, orig.ident))
write.table(meta.data, "./data/raw/bird/metadata.csv", sep=",")
write.table(counts, "./data/raw/bird/counts.csv", sep=",")

# Create Seurat
# "Row names in the metadata need to match the column names of the counts matrix"
# Metadata cells x vars. Counts: genes x cells
counts <- t(counts) # cells x genes -> genes x cells
colnames(counts) <- meta.data$cell
rownames(meta.data) <- meta.data$cell
meta.data <-meta.data %>% select(!cell)
bird <- CreateSeuratObject(counts = counts, project = "bird", meta.data = meta.data)
bird # 19857 features across 29497 samples within 1 assay 
SaveH5Seurat(bird, "./data/raw/bird/bird")


