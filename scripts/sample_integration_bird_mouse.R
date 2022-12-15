# # https://satijalab.org/seurat/articles/integration_introduction.html
library(dplyr)
library(Seurat)
library(SeuratDisk)

# Integrated analysis: standard workflow 
integrated.analysis <- function(combined.seurat){
  combined.seurat <- ScaleData(combined.seurat, verbose = TRUE)
  combined.seurat <- RunPCA(combined.seurat, npcs = 30, verbose = TRUE)
  combined.seurat <- RunUMAP(combined.seurat, reduction = "pca", dims = 1:30)
  combined.seurat <- FindNeighbors(combined.seurat, reduction = "pca", dims = 1:30)
  combined.seurat <- FindClusters(combined.seurat, resolution = 0.5)
  return(combined.seurat)
}

for (selected_class in c('GABAergic', 'Glutamatergic')){
  setwd("/mnt/data/joram/elfn1_evolution")
  
  print(paste("Processing", selected_class))
  # Load SCTransformed data
  mouse <- LoadH5Seurat("./data/raw/mouse/mouseSCT.h5seurat") 
  mouse <- mouse[,mouse@meta.data$class == selected_class]
  gc()
  # Next, bird
  bird <- LoadH5Seurat("./data/raw/bird/birdSCT.h5seurat")
  bird@meta.data <- bird@meta.data %>% 
    rename(brain_region = position2, organism = species, cluster = cluster_int_sub2)
  # Add class (GABAergic, Glutamatergic etc)
  class <- bird@meta.data$cluster
  class[startsWith(class, "GABA")] <- "GABAergic"
  class[startsWith(class, "RA_Glut")] <- "Glutamatergic"
  class[startsWith(class, "HVC_Glut")] <- "Glutamatergic"
  #class[startsWith(class, "Pre")] <- "Glutamatergic"
  bird@meta.data['class'] <- class
  # Select neurons
  bird <- bird[,bird@meta.data$class == selected_class]
  
  # split across areas and species. Only need brain_region for glut 
  object.list <- SplitObject(bird, split.by = "brain_region")
  rm(bird)
  gc()
  object.list['mouse'] <- mouse
  rm(mouse)
  
  # Integrate: find anchors. Exclude first if using Glu
  if (selected_class == "Glutamatergic"){
    # delete 18 samples from RA - too small for integration
    object_list <- object_list[-1]
  }
  anchors <- FindIntegrationAnchors(object.list = object.list[-1])
  #rm(object.list)
  # Do the actual integration
  object.combined <- IntegrateData(anchorset = anchors)
  gc()
  
  # Usual integrated analysis
  object.combined <- integrated.analysis(object.combined)
  
  # Visualization
  DimPlot(object.combined, reduction = "umap", group.by = "organism")
  
  fname <- paste0("mouse_bird_", selected_class, "_integrated")
  savename <- paste0("./data/raw/", fname)
  SaveH5Seurat(object.combined, savename, overwrite=TRUE)
  setwd("./data/raw/")
  # Also save as H5ad for viz in scanpy
  Convert(paste0(fname, ".h5seurat"), dest = "h5ad")
  
  # Plot
  p1 <- DimPlot(object.combined, reduction = "umap", group.by = "subclass")
  p2 <- DimPlot(object.combined, reduction = "umap", group.by="organism")
  p1 + p2
  # clear memory
  rm(list=ls())
  gc()
}