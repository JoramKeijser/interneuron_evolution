# # https://satijalab.org/seurat/articles/integration_introduction.html
library(dplyr)
library(Seurat)
library(SeuratDisk)

bird <- LoadH5Seurat("./data/raw/bird/birdSCT.h5seurat")


# Split by area
bird.list <- SplitObject(bird, split.by = "position2")
# 1. Integrate each area
#https://satijalab.org/seurat/reference/FindIntegrationAnchors.html
bird.anchors <- FindIntegrationAnchors(object.list = bird.list)
# this command creates an 'integrated' data assay
bird.combined <- IntegrateData(anchorset = bird.anchors)
# 2. Integrate area and species
# split across areas and species
bird.list <- SplitObject(bird.combined, split.by = "position2")
bird.list <- sapply(bird.list, FUN = function(x){SplitObject(x, split.by="species")})
# Find new anchors, 
bird.anchors <- FindIntegrationAnchors(object.list = bird.list)
bird.combined.twice <- IntegrateData(anchorset = bird.anchors)

# Now perform integrated analysis
# Run the standard workflow for visualization and clustering
integrated.analysis <- function(combined.seurat){
  combined.seurat <- ScaleData(combined.seurat, verbose = TRUE)
  combined.seurat <- RunPCA(combined.seurat, npcs = 30, verbose = TRUE)
  combined.seurat <- RunUMAP(combined.seurat, reduction = "pca", dims = 1:30)
  combined.seurat <- FindNeighbors(combined.seurat, reduction = "pca", dims = 1:30)
  combined.seurat <- FindClusters(combined.seurat, resolution = 0.5)
  return(combined.seurat)
}

bird.combined <- integrated.analysis(bird.combined)
Idents(object = bird.combined) <- bird.combined@meta.data$cluster_int
# Visualization
p1 <- DimPlot(bird.combined, reduction = "umap", group.by = "position2")
p2 <- DimPlot(bird.combined, reduction = "umap", group.by="species")
p1 + p2

# do the same but for twice combined
bird.combined.twice <- integrated.analysis(bird.combined.twice)
Idents(object = bird.combined.twice) <- bird.combined.twice@meta.data$cluster_int
p1 <- DimPlot(bird.combined.twice, reduction = "umap", group.by = "position2")
p2 <- DimPlot(bird.combined.twice, reduction = "umap", group.by="species")
p3 <- DimPlot(bird.combined.twice, reduction = "umap", group.by="ident")
p1 + p2 + p3


#mouse <- LoadH5Seurat("./data/raw/mouse/mouseSCT.h5seurat")
#colnames(mouse@meta.data)
