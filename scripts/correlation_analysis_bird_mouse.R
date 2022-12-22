# Correlate gene expression between clusters
library(dplyr)
library(Seurat)
library(SeuratDisk)
library(gplots)

source("./src/seurat_utils.R")

cell_type <- "Gaba" # "Gaba" or "Glut

# Load bird data
bird <- LoadH5Seurat("./data/Colquitt21/bird.h5seurat") 
bird # 218342 features across 29497 samples within 1 assay 
# Select largest dataset
bird <- bird[, bird@meta.data$species == 'zf']
# 18342 features across 25392 samples within 1 assay 
# Mouse data
mouse <- LoadH5Seurat("./data/Tasic18/mouse.h5seurat")
mouse # 42094 features across 15413 samples within 1 assay 

# 0. Select only GABAergic or glutametergic cells
if (cell_type == "Gaba"){
  bird <- bird[, bird@meta.data$cluster_int_sub2 %>% startsWith("GABA")]  # 4956 samples
  mouse <- mouse[, mouse@meta.data$class == "GABAergic"] # 6125 samples
} else if (cell_type == "Glut"){
  turtle <- turtle[, turtle@meta.data$clusters %>% startsWith("e")]  # 640 samples
  mouse <- mouse[, mouse@meta.data$class == "Glutamatergic"] # 7366 samples
}
# Merge bird clusters
bird@meta.data <- bird@meta.data %>% mutate(subclass = recode(cluster_int_sub2, 
                          "GABA-1-1" = "1", "GABA-1-2" = "1",
                          "GABA-2" = "2", "GABA-3" = "3", "GABA-4" = "4", 
                          "GABA-5-1" = "5", "GABA-5-2" = "5", "GABA-5-3" = "5",
                          "GABA-6" = "6", "GABA-7" = "7", "GABA-8"="8",
                          "GABA-Pre" = "Pre", ))

#mouse <- clustering_pipeline(mouse, clustering_resolution = 0.4) 
#bird <- clustering_pipeline(bird, clustering_resolution=0.4) 
# Which groups/clusters to compare?
Idents(bird) <- bird@meta.data$subclass
Idents(mouse) <- mouse@meta.data$subclass

# 1. For each dataset: marker genes findAllMarkers()
# To do: courser clustering of turtle neurons? Now have 17. 
# Or finer clustering of mouse neurons. 
top.bird.markers <- get_top_markers(bird) # 2329, 2343
top.mouse.markers <- get_top_markers(mouse) # Gaba: 2588, Glut: 4798

# 3. Intersect to select shared marker genes
shared.markers <- intersect(top.mouse.markers, top.bird.markers)
print(length(shared.markers)) # ~474

# 4. Average across clusters. slot = 'counts' -> don't exponentiate
# 5. Log(x+1) + 0.1
# 6. : across-cluster mean
bird.avg <- score_genes(bird, top_markers = shared.markers)
mouse.avg <- score_genes(mouse, top_markers = shared.markers)

# 6. Correlate. 
C <- cor(bird.avg, mouse.avg) # 18 x 7
heatmap(C)
# Bird - delete 7,8,Pre. 
# Mouse - delete Serpinf1, Sncg
if (cell_type == "Gaba"){
  C <- as.data.frame(C)
  Csub <- C[!(rownames(C) %in% c("7", "8", "Pre")), !(colnames(C) %in% c("Serpinf1", "Sncg"))]
} else if (cell_type == "Glut"){
  # to do: subset 
}

# Fixed color range
colors <- seq(-.25,.35,length=100)
colormap <- hcl.colors(99, "RdYlBu",rev=TRUE)
heatmap.2(as.matrix(Csub), col = colormap, 
         breaks=colors, density.info="none", trace="none", 
         dendrogram=c("row"), symm=F,symkey=F,symbreaks=T)

# 7. Redo after shuffling to get significance. 


