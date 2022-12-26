# Standard workflow on human data (Bakken et al)
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)

human <- LoadH5Seurat("./data/Bakken21/bakken.h5seurat")

# QC
mito.count <- colSums(human[sapply(rownames(human), function(x) {startsWith(x, "MT") }), ])
total.count <- colSums(human) # == nCount_RNA
human[['percent.mito']] <- mito.count / total.count * 100
VlnPlot(human, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
# percent.mito is super low

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(human, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(human, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# to do: subset
human <- NormalizeData(object = human) # log norm
human <- FindVariableFeatures(object = human) # HVGs 
human <- ScaleData(human) #
human <- RunPCA(object = human) 
human <- FindNeighbors(object = human)
human <- FindClusters(object = human, resolution = 0.4)
human <- RunUMAP(object = human, dims=1:30)

# Viz. 
DimPlot(object = human, reduction = "umap")
FeaturePlot(human, c("SST", "PVALB", "VIP", "LAMP5"))

# PCA plot
DimPlot(object = human, reduction = "pca")
FeaturePlot(human, c("SST", "PVALB", "VIP", "LAMP5"), reduction = "pca")
