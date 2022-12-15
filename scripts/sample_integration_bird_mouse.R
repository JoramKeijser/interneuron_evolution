# # https://satijalab.org/seurat/articles/integration_introduction.html
library(dplyr)
library(Seurat)
library(SeuratDisk)

bird <- LoadH5Seurat("./data/raw/bird/birdSCT.h5seurat")

bird@meta.data <- bird@meta.data %>% 
  rename(brain_region = position2, organism = species, cluster = cluster_int_sub2)
# Add class GABAergic/...
subclass <- bird@meta.data$cluster
subclass[startsWith(subclass, "GABA")] <- "GABA"
subclass[startsWith(subclass, "RA_Glut")] <- "Glut"
subclass[startsWith(subclass, "HVC_Glut")] <- "Glut"
subclass[startsWith(subclass, "Pre")] <- "Pre"
subclass

class <- subclass
class[class == "GABA"] <- "neuron"
class[class == "Glut"] <- "neuron"
#bird@meta.data['class'] <- "GABAergic"
bird@meta.data$class <- class
# Select
#bird <- bird[,bird@meta.data['subclass']  == "GABA"]
# OR:
bird <- bird[,subclass == "Glut"]
bird@meta.data['class']
# split across areas and species
object.list <- SplitObject(bird, split.by = "brain_region")
object.list <- sapply(object.list, FUN = function(x){SplitObject(x, split.by="organism")})
object.list
rm(bird)

# Load mouse
mouse <- LoadH5Seurat("./data/raw/mouse/mouseSCT.h5seurat") 
# check gene names
length(intersect(rownames(mouse), rownames(object.list[[1]]))) # 12040
#mouse <- mouse[,mouse@meta.data$class == "GABAergic"]
# Or: 
mouse <- mouse[,mouse@meta.data$class == "Glutamatergic"]
object.list['mouse'] <- mouse
rm(mouse)

# Integrate: find anchors. Exclude first if using Glu
anchors <- FindIntegrationAnchors(object.list = object.list[-1])
# Do the actual integration
object.combined <- IntegrateData(anchorset = anchors)

# Integrated analysis: standard workflow 
integrated.analysis <- function(combined.seurat){
  combined.seurat <- ScaleData(combined.seurat, verbose = TRUE)
  combined.seurat <- RunPCA(combined.seurat, npcs = 30, verbose = TRUE)
  combined.seurat <- RunUMAP(combined.seurat, reduction = "pca", dims = 1:30)
  combined.seurat <- FindNeighbors(combined.seurat, reduction = "pca", dims = 1:30)
  combined.seurat <- FindClusters(combined.seurat, resolution = 0.5)
  return(combined.seurat)
}

object.combined <- integrated.analysis(object.combined)
#Idents(object = object.combined) <- object.combined@meta.data$cluster_int
# Visualization
DimPlot(object.combined, reduction = "umap", group.by = "organism")
SaveH5Seurat(object.combined, "./data/raw/mouse_bird_glut_integrated")
setwd("./data/raw/")
Convert("mouse_bird_glut_integrated.h5seurat", dest = "h5ad")
gaba.combined <- LoadH5Seurat("mouse_bird_gaba_integrated.h5seurat") 
Convert("mouse_bird_gaba_integrated.h5seurat", dest = "h5ad")

p1 <- DimPlot(object.combined, reduction = "umap", group.by = "")
p2 <- DimPlot(object.combined, reduction = "umap", group.by="species")
p1 + p2

# Subset GABA and Gluta cells. 