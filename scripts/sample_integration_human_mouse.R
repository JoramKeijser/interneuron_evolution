# # https://satijalab.org/seurat/articles/integration_introduction.html
library(dplyr)
library(Seurat)
library(SeuratDisk)
source("./src/seurat_utils.R") # integrated.analysis pipeline
# To do: use SCTransform or not?
use_SCT <- FALSE
datadir <- "./data/"

#setwd("/mnt/data/joram/elfn1_evolution")
  
# Load data
mouse <- LoadH5Seurat(paste0(datadir, "Tasic18/", "mouse.h5seurat"))
mouse <- mouse[,mouse@meta.data$class == "GABAergic"]
# Same naming convention
mouse@meta.data <- mouse@meta.data %>% mutate(percent.mito = percent_mt_exon_reads)
# Next, human
human <- LoadH5Seurat("./data/Bakken21/bakken.h5seurat")
mito.count <- colSums(human[sapply(rownames(human), function(x) {startsWith(x, "MT") }), ])
total.count <- colSums(human) # == nCount_RNA
human[['percent.mito']] <- mito.count / total.count * 100
human[['organism']] <- 'homo sapiens'
object.list = list(human = human, mouse = mouse)
if (use_SCT){
  # SCT workflow
  object.list <- lapply(object.list, 
                        function(x){SCTransform(x, vars.to.regress = "percent.mito")})
  features <- SelectIntegrationFeatures(object.list, nfeatures=3000)
  object.list <- PrepSCTIntegration(object.list, anchor.features = features)
  anchors <- FindIntegrationAnchors(object.list, 
                                    normalization.method = 'SCT', 
                                    anchor.features = features)
  combined <- IntegrateData(anchors, normalization.method = 'SCT')
  combined <- Integrated.analysis(combined, FALSE)  
} else{
  object.list <- lapply(object.list, function(x){ 
                            x %>% NormalizeData() %>% FindVariableFeatures()})
  anchors <- FindIntegrationAnchors(object.list = object.list)
  combined <- IntegrateData(anchorset = anchors)
  # Usual integrated analysis
  combined <- integrated.analysis(combined, do_scale = TRUE)
}

  
# Visualization
DimPlot(combined, reduction = "umap", group.by = "organism")
# Together with cell types from mouse
p1 <- DimPlot(combined, reduction = "umap", group.by = "subclass")
p2 <- DimPlot(combined, reduction = "umap", group.by="organism")
p1 + p2

fname <- paste0("mouse_human_GABAergic", "_integrated")
getwd()
savename <- paste0("./data/integrated_datasets/", fname)
SaveH5Seurat(combined, savename, overwrite=TRUE)
# Also save as H5ad for viz in scanpy
Convert(paste0(savename, ".h5seurat"), dest = "h5ad", overwrite=TRUE)
