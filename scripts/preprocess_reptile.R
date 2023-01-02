# Save raw bird data as Seurat
library(dplyr)
library(Seurat)
library(SeuratDisk)
rm(list = ls())
gc()

# Where to load from and save to
projectdir <- "/home/joram/Dropbox/elfn1_evolution"
if (getwd() != projectdir){
  setwd(projectdir)  
}
datadir <- "./data/raw/tosches/"
savedir <- datadir
species_list <- c("turtle", "lizard")

for (species in species_list){
  # Load data    
  load(paste0(datadir, paste0(species, ".neurons.Robj")), verbose=TRUE)
  
  # use same name ("neurons") for both species, to avoid further if-statements
  if (species == "turtle"){
    neurons <- turtle.neurons
    rm(turtle.neurons)
  } else if (species == "lizard"){
    neurons <- lizard.neurons
    rm(lizard.neurons)
  }
  
  # Extract data
  counts <- attr(neurons, "raw.data") # genes x cells 
  meta.data <- attr(neurons, "data.info") # cells x features 
  
  # Gene names in "Elfn1" format (mouse convention)
  rownames(counts) <- sapply(rownames(counts), function(name){
    paste0(substr(name, 1, 1), tolower(substr(name, 2, nchar(name))))
  })

  # Create and save Seurat object
  neurons <- CreateSeuratObject(counts = counts, project = species,
                               meta.data = meta.data, min.features = 1)
  print(neurons) # 23500 features across 8575 samples within 1 assay
  print(paste0("Save Seurat object to ", savedir, species))
  SaveH5Seurat(neurons, paste0(savedir, species), overwrite=TRUE)
  
  # Also save in H5ad/anndata format
  Convert(paste0(savedir, species, ".h5seurat"), 
          dest = "h5ad", overwrite=TRUE)
}



