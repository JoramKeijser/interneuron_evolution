# Save raw bird data as Seurat

library(dplyr)
library(Seurat)
library(SeuratDisk)
rm(list = ls())
gc()

# Where to load from and save to
datadir <- "~/Dropbox/scRNAseq_data/Tosches18/"

for (species in c("turtle", "lizard")){
  print(paste("Processing", species))
  savedir <- paste0("./data/Tosches18_", species, "/")
  # Load data    
  load(paste0(datadir, paste0(species, ".neurons.Robj")), verbose=TRUE)
  
  # give it the same name
  if (species == "turtle"){
    neurons <- turtle.neurons
    rm(turtle.neurons)
  } else if (species == "lizard"){
    neurons <- lizard.neurons
    rm(lizard.neurons)
  }
  
  # Extract data
  data.attributes <- attributes(neurons)
 
  counts <- attr(neurons, "raw.data") # genes x cells
  cell.types <- attr(neurons, "ident")
  cell.names <- attr(neurons, "cell.names")
  gene.names <- rownames("counts")
  data.info <- attr(neurons, "data.info")
  
  # Save it
  print(paste("Save to ", savedir))
  write.table(data.info, paste0(savedir, "metadata.csv"), sep=",")
  write.table(counts, paste0(savedir, "counts.csv"), sep=",")
  
  # Create Seurat
  seurat_object <- CreateSeuratObject(counts = counts, project = species,
                               meta.data = data.info, 
                               min.cells = 1, min.features = 1)
  SaveH5Seurat(seurat_object, paste0(savedir, species))
  # Also save in H5ad format for scanpy
  Convert(paste0(savedir, paste0("mouse", ".h5seurat")), 
          dest = "h5ad", overwrite=TRUE)
}



