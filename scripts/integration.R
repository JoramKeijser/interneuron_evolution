# Integrate bird and mouse samples
projectdir <- "/home/joram/Dropbox/elfn1_evolution" # update this
setwd(projectdir)  
# Directories where to load and save data 
datadir <- "./data/seurat/"
savedir_seurat <- "./data/seurat/"
savedir_anndata <- "./data/anndata/"
source("./src/seurat_utils.R") # integrated.analysis & packages
use_SCT <- TRUE 

# Separately treat exc and inh neurons
for (selected_class in c('GABAergic',  'Glutamatergic')){
  print(paste("Processing", selected_class))
  # Load data
  mouse <- LoadH5Seurat(paste0(datadir, "mouse.h5seurat"))
  mouse <- mouse[,mouse@meta.data$class == selected_class]
  # Same naming convention
  mouse@meta.data <- mouse@meta.data %>% mutate(percent.mito = percent_mt_exon_reads)
  # Next, bird
  bird <- LoadH5Seurat(paste0(datadir,"bird.h5seurat"))
  # Only use largest dataset (zebra finch (zf), not Bengalese finch (bf))
  bird <- bird[, bird@meta.data$species == 'zf'] # 25392 samples within 1 assay 
  bird@meta.data <- bird@meta.data %>% 
    rename(brain_region = position2, organism = species, cluster = cluster_int_sub2)
  # Select neuron class
  if (selected_class == "GABAergic"){
    bird <- bird[,sapply(bird@meta.data$cluster, function(x) {grepl("GABA", x)})]
  } else if (selected_class == "Glutamatergic"){
    bird <- bird[,sapply(bird@meta.data$cluster, function(x) {grepl("Glut", x)})]
  }
  
  # split across areas and species. Only need brain_region for glut 
  object.list <- c("songbird" = bird, "mouse" = mouse) #SplitObject(bird, split.by = "organism")
  # Only 1 mouse region
  #object.list['mouse'] <- mouse
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
    combined <- integrated.analysis(combined, scale = FALSE)  
    fname <- paste0("mouse_bird_", selected_class, "_integrated_SCT")
  } else{
    object.list <- lapply(object.list, function(x){ 
                              x %>% NormalizeData() %>% FindVariableFeatures()})
    anchors <- FindIntegrationAnchors(object.list = object.list)
    combined <- IntegrateData(anchorset = anchors)
    # Usual integrated analysis
    combined <- integrated.analysis(combined, scale = TRUE)
    fname <- paste0("mouse_bird_", selected_class, "_integrated")
  }

  # Visualization
  p1 <- DimPlot(combined, reduction = "umap", group.by = "subclass")
  p2 <- DimPlot(combined, reduction = "umap", group.by="organism")
  p1 + p2
  
  # Save Seurat object
  savename <- paste0(savedir_seurat, fname)
  SaveH5Seurat(combined, savename, overwrite=TRUE)
  # Also save as H5ad for viz in scanpy
  Convert(paste0(savedir_seurat, paste0(fname, ".h5seurat")), 
          dest = paste0(savedir_anndata, fname, ".h5ad"), overwrite=TRUE)
  print("")
}