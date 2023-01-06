# Correlate gene expression between clusters
source("./src/seurat_utils.R")

for (cell_type in c("Glut", "GABA")){
  # Load bird data
  bird <- LoadH5Seurat("./data/seurat/bird.h5seurat") 
  bird # 218342 features across 29497 samples within 1 assay 
  # Select largest dataset: zebra finch (not bengalese)
  bird <- bird[, bird@meta.data$species == 'zf'] # 25392 samples within 1 assay 
  # Load mouse data
  mouse <- LoadH5Seurat("./data/seurat/tasic/mouse.h5seurat")
  mouse # 42094 features across 15413 samples within 1 assay 
  
  # 0. Select only GABAergic or glutametergic cells
  if (cell_type == "GABA"){
    bird <- bird[, bird@meta.data$cluster_int_sub2 %>% startsWith("GABA")]  # 640 samples
    mouse <- mouse[, mouse@meta.data$class == "GABAergic"] # 6125 samples
    # Merge bird clusters into coarser subclasses 
    bird@meta.data <- bird@meta.data %>% 
      mutate(subclass = recode(cluster_int_sub2, 
                      "GABA-1-1" = "1", "GABA-1-2" = "1",
                       "GABA-2" = "2", "GABA-3" = "3", "GABA-4" = "4", 
                       "GABA-5-1" = "5", "GABA-5-2" = "5", "GABA-5-3" = "5",
                       "GABA-6" = "6", "GABA-7" = "7", "GABA-8"="8",
                       "GABA-Pre" = "Pre", ))

  } else if (cell_type == "Glut"){
    bird <- bird[,sapply(bird@meta.data$cluster_int_sub2, function(x){grepl("Glut", x)})]  # 4956 samples
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
  print(paste(length(shared.markers), "shared marker genes")) 
  
  # Average across clusters; log1p; across-cluster mean
  bird.avg <- score_genes(bird, top_markers = shared.markers)
  mouse.avg <- score_genes(mouse, top_markers = shared.markers)
  
  # Correlate
  C <- cor(mouse.avg, bird.avg) 
  # Throw out un-matched songbird cells. 
  if (cell_type == "GABA"){
    C <- C[!rownames(C) %in% c("Serpinf1", "Sncg"), !colnames(C) %in% c("Pre", "7", "8")]
    # Reorder
    C <- C[c("Meis2", "Sst", "Pvalb", "Vip", "Lamp5"),]
    }
  else if (cell_type == "Glut"){
    C <- C[, !(colnames(C) %in% c("NP", "CR"))]
    colnames(C) <- sapply(colnames(C), function(x){str_replace(x, "_Glut", "")})
  }
  # Plot it
  fontsize <- 25
  # consistent color range for GABA and Glut:
  mat_breaks <- seq(-0.3, 0.4, length.out = 100) 
  pdf(file = paste0("./figures/corr_", cell_type, ".pdf")) # save
  setHook("grid.newpage", 
          function() pushViewport(viewport(x=1,y=1, width=.9, height=0.7, 
                                           name="vp", just=c("right","top"))), action="prepend")
  pheatmap(C,  legend=TRUE, main = "RNA correlation", breaks=mat_breaks,
           legend_breaks = c(-.3, 0, 0.3),  legend_labels = c("-0.3", "0", "0.3"),
           fontsize=fontsize, srtCol=45)
  setHook("grid.newpage", NULL, "replace")
  grid.text("mouse clusters", x=-0.07, y = 0.4, rot=90, gp=gpar(fontsize=fontsize))
  grid.text("songbird clusters", y=-0.07, x=.4, gp=gpar(fontsize=fontsize))
  dev.off()
}



