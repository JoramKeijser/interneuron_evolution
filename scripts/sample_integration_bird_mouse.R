# # https://satijalab.org/seurat/articles/integration_introduction.html
library(dplyr)
library(Seurat)
library(SeuratDisk)
source("./src/seurat_utils.R") # integrated.analysis pipeline
# To do: use SCTransform or not?
use_SCT <- FALSE
datadir <- "./data/"


for (selected_class in c('GABAergic', 'Glutamatergic')){
  setwd("/mnt/data/joram/elfn1_evolution")
  
  print(paste("Processing", selected_class))
  # Load data
  mouse <- LoadH5Seurat(paste0(datadir, "Tasic18/", "mouse.h5seurat"))
  mouse <- mouse[,mouse@meta.data$class == selected_class]
  # Same naming convention
  mouse@meta.data <- mouse@meta.data %>% mutate(percent.mito = percent_mt_exon_reads)
  # Next, bird
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
  bird <- bird[,bird@meta.data$class == selected_class]
  
  # split across areas and species. Only need brain_region for glut 
  object.list <- SplitObject(bird, split.by = "brain_region")
  if (selected_class == "Glutamatergic"){
    # delete 18 samples from RA - too small for integration
    object_list <- object_list[-1]
  }
  object.list['mouse'] <- mouse
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
    combined <- integrated.analysis(combined, scale = TRUE)
  }

  # Visualization
  DimPlot(combined, reduction = "umap", group.by = "organism")
  # Together with cell types from mouse
  p1 <- DimPlot(combined, reduction = "umap", group.by = "subclass")
  p2 <- DimPlot(combined, reduction = "umap", group.by="organism")
  p1 + p2
  
  fname <- paste0("mouse_bird_", selected_class, "_integrated")
  savename <- paste0("./data/raw/", fname)
  SaveH5Seurat(combined, savename, overwrite=TRUE)
  setwd("./data/raw/")
  # Also save as H5ad for viz in scanpy
  Convert(paste0(fname, ".h5seurat"), dest = "h5ad")
  
 
  # clear memory
  gc()
}