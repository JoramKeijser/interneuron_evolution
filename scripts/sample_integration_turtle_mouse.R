# # https://satijalab.org/seurat/articles/integration_introduction.html
library(dplyr)
library(Seurat)
library(SeuratDisk)
source(seurat_utils) # integrated.analysis pipeline
# To do: use SCTransform or not?
use_SCT <- FALSE
datadir <- "./data/"

for (selected_class in c('GABAergic', 'Glutamatergic')){
  
  print(paste("Processing", selected_class))
  # Load data
  mouse <- LoadH5Seurat(paste0(datadir, "Tasic18/", "mouse.h5seurat"))
  mouse <- mouse[,mouse@meta.data$class == selected_class]
  # Same naming convention
  mouse@meta.data <- mouse@meta.data %>% mutate(percent.mito = percent_mt_exon_reads)
  # Next, bird
  turtle <- LoadH5Seurat(paste0(datadir,"Tosches18_turtle/",  "turtle.h5seurat"))
  turtle@meta.data <- turtle@meta.data %>% 
    mutate(brain_region = areaident, 
           percent.mito = pMITO / 100, 
           cluster = clusters)
  turtle@meta.data$organism <- "turtle" 
  # Add class (GABAergic, Glutamatergic etc
  if (selected_class == "GABAergic"){
    turtle <- turtle[, turtle@meta.data$cluster %>% startsWith("i")]
  } else if (selected_class == "Glutamatergic"){
    turtle <- turtle[, turtle@meta.data$cluster %>% startsWith("e")]
  }

  object.list <- list(mouse = mouse, turtle = turtle)
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
  
  
  # Save integrated object
  if (use_SCT){
    fname <- paste0("mouse_turtle_", selected_class, "_integrated_SCT")  
  } else {
    fname <- paste0("mouse_turtle_", selected_class, "_integrated")  
  }
  savename <- paste0("./data/integrated_datasets/", fname)
  SaveH5Seurat(combined, savename, overwrite=TRUE)
  # Also save as H5ad for viz in scanpy
  Convert(paste0(savename, ".h5seurat"), dest = "h5ad", overwrite=TRUE)
  
  # Plot
  p1 <- DimPlot(combined, reduction = "umap", group.by = "subclass")
  p2 <- DimPlot(combined, reduction = "umap", group.by="organism")
  p1 + p2
}