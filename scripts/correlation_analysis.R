# Correlate gene expression between clusters
library(dplyr)
library(Seurat)
library(SeuratDisk)
library(gplots)
library(stringr)
rm(list = ls())

projectdir <- "/home/joram/Dropbox/elfn1_evolution"
setwd(projectdir)  
source("./src/seurat_utils.R")


for (cell_type in c("Glut", "GABA")){
  # Load bird data
  bird <- LoadH5Seurat("./data/seurat/bird.h5seurat") 
  bird # 218342 features across 29497 samples within 1 assay 
  # Select largest dataset: zebra finch (not bengalese)
  bird <- bird[, bird@meta.data$species == 'zf'] # 25392 samples within 1 assay 
  # Load mouse data
  mouse <- LoadH5Seurat("./data/raw/tasic/mouse.h5seurat")
  mouse # 42094 features across 15413 samples within 1 assay 
  
  # 0. Select only GABAergic or glutametergic cells
  if (cell_type == "GABA"){
    bird <- bird[, bird@meta.data$cluster_int_sub2 %>% startsWith("GABA")]  # 4956 samples
    mouse <- mouse[, mouse@meta.data$class == "GABAergic"] # 6125 samples
    # Merge bird clusters into coarser subclasses 
    bird@meta.data <- bird@meta.data %>% mutate(subclass = 
                                                  recode(cluster_int_sub2, 
                                                         "GABA-1-1" = "1", "GABA-1-2" = "1",
                                                         "GABA-2" = "2", "GABA-3" = "3", "GABA-4" = "4", 
                                                         "GABA-5-1" = "5", "GABA-5-2" = "5", "GABA-5-3" = "5",
                                                         "GABA-6" = "6", "GABA-7" = "7", "GABA-8"="8",
                                                         "GABA-Pre" = "Pre", ))

  } else if (cell_type == "Glut"){
    bird <- bird[,sapply(bird@meta.data$cluster_int_sub2, function(x){grepl("Glut", x)})]  # 640 samples
    bird@meta.data  <- bird@meta.data %>% mutate(subclass = cluster_int_sub2)
    mouse <- mouse[, mouse@meta.data$class == "Glutamatergic"] # 7366 samples
  }
  # Which groups/clusters to compare?
  Idents(bird) <- bird@meta.data$subclass
  Idents(mouse) <- mouse@meta.data$subclass
  
  # For each dataset: marker genes 
  top.bird.markers <- get_top_markers(bird) 
  top.mouse.markers <- get_top_markers(mouse) 
  
  # Select shared marker genes
  shared.markers <- intersect(top.mouse.markers, top.bird.markers)
  print(paste(length(shared.markers), "shared marker genes")) # 475-500
  
  # Average across clusters; log1p; across-cluster mean
  bird.avg <- score_genes(bird, top_markers = shared.markers)
  mouse.avg <- score_genes(mouse, top_markers = shared.markers)
  
  # 6. Correlate. 
  C <- cor(bird.avg, mouse.avg) # 18 x 7
  C <- as.data.frame(C)
  # Bird - delete 7,8,Pre. 
  # Mouse - delete Serpinf1, Sncg
  if (cell_type == "Gaba"){
    # Remove small subclasses or those that correspond to subcortical INs
    C <- C[!(rownames(C) %in% c("7", "8", "Pre")), !(colnames(C) %in% c("Serpinf1", "Sncg"))]
    # Reorder
  } else if (cell_type == "Glut"){
    C <- C[, !(colnames(C) %in% c("NP", "CR"))]
    rownames(C) <- sapply(rownames(C), function(x){str_replace(x, "_Glut", "")})
  }
  
  # Fixed color range
  colors <- seq(-.3,.35,length=100)
  colormap <- hcl.colors(99, "RdYlBu",rev=TRUE)
  # Make and save the colormap
  pdf(file = paste0("./figures/corr_", cell_type, ".pdf"))
  heatmap.2(as.matrix(C), col = colormap, 
            breaks=colors, density.info="none", trace="none", 
            dendrogram=c("both"), symm=F,symkey=F,symbreaks=T, 
            srtCol=45, srtRow = 0, key.title = "correlation", key.xlab = "")
  
  dev.off()
  # Save
  write.table(C, paste0("./data/correlation_", cell_type, ".csv"))
  write.table(bird.avg, paste0("./data/marker_expr_bird_", cell_type, ".csv"))
  write.table(mouse.avg, paste0("./data/marker_expr_mouse_ ", cell_type, ".csv"))
  # Todo: shuffling to get significance. 
}



