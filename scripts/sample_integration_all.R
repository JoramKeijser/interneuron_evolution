# # https://satijalab.org/seurat/articles/integration_introduction.html
library(dplyr)
library(Seurat)
library(SeuratDisk)
source("./src/seurat_utils.R") # integrated.analysis pipeline
# To do: use SCTransform or not?
use_SCT <- TRUE
datadir <- "./data/"

#setwd("/mnt/data/joram/elfn1_evolution")
  
# Load data
mouse <- LoadH5Seurat(paste0(datadir, "Tasic18/", "mouse.h5seurat"))
mouse <- mouse[,mouse@meta.data$class == "GABAergic"]
# Same naming convention
mouse@meta.data <- mouse@meta.data %>% mutate(percent.mito = percent_mt_exon_reads)
# Next, human
source("./scripts/human2seurat.R")
# Bird
datadir <- "./data/"
bird <- LoadH5Seurat(paste0(datadir,"Colquitt21/",  "bird.h5seurat"))
bird@meta.data <- bird@meta.data %>% 
  rename(brain_region = position2, organism = species, cluster = cluster_int_sub2)
# Add class (GABAergic, Glutamatergic etc)
class <- bird@meta.data$cluster
class[startsWith(class, "GABA")] <- "GABAergic"
class[startsWith(class, "RA_Glut")] <- "Glutamatergic"
class[startsWith(class, "HVC_Glut")] <- "Glutamatergic"
bird@meta.data['class'] <- class
# Select neurons
bird <- bird[,bird@meta.data$class == "GABAergic"]
bird <- bird[, bird@meta.data$organism == "zf"]
# Finally, turtle
turtle <- LoadH5Seurat(paste0(datadir,"Tosches18_turtle/",  "turtle.h5seurat"))
turtle@meta.data <- turtle@meta.data %>% 
  mutate(brain_region = areaident, 
         percent.mito = pMITO / 100, 
         cluster = clusters)
turtle@meta.data[['organism']] <- "turtle" 
# Add class (GABAergic, Glutamatergic etc
turtle <- turtle[, turtle@meta.data$cluster %>% startsWith("i")]

# Combine all datasets
object.list = list(human = human, mouse = mouse, bird = bird, turtle = turtle)
gc()
# Normalization, HVGs
object.list <- lapply(X = object.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

# Select features, do PCA
features <- SelectIntegrationFeatures(object.list = object.list)
object.list <- lapply(X = object.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# Integrate, with mouse data as reference
anchors <- FindIntegrationAnchors(object.list = object.list, reference = c(2), reduction = "rpca",
                                  dims = 1:30)
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

# Standard work flow
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, verbose = FALSE)
integrated <- RunUMAP(integrated, dims = 1:30)

# Visualization
DimPlot(integrated, reduction = "umap", group.by = "organism")
# Together with cell types from mouse
p1 <- DimPlot(integrated, reduction = "umap", group.by = "class")
p2 <- DimPlot(integrated, reduction = "umap", group.by="organism")
p1 + p2


fname <- paste0("GABAergic_integrated")
savename <- paste0("./data/integrated_datasets/", fname)
SaveH5Seurat(integrated, savename, overwrite=TRUE)
# Also save as H5ad for viz in scanpy
Convert(paste0(savename, ".h5seurat"), dest = "h5ad", overwrite=TRUE)

