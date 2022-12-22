# Correlate gene expression between clusters
library(dplyr)
library(Seurat)
library(SeuratDisk)
library(gplots)

source("./src/seurat_utils.R")

cell_type <- "Gaba" # "Gaba" or "Glut

# Load turtle data
turtle <- LoadH5Seurat("./data/Tosches18_turtle/turtle.h5seurat") 
turtle # 22915 features across 8575 samples within 1 assay 
mouse <- LoadH5Seurat("./data/Tasic18/mouse.h5seurat")
mouse # 42094 features across 15413 samples within 1 assay 
head(mouse@meta.data)

# 0. Select only GABAergic or glutametergic cells
if (cell_type == "Gaba"){
  turtle <- turtle[, turtle@meta.data$clusters %>% startsWith("i")]  # 4956 samples
  mouse <- mouse[, mouse@meta.data$class == "GABAergic"] # 6125 samples
} else if (cell_type == "Glut"){
  turtle <- turtle[, turtle@meta.data$clusters %>% startsWith("e")]  # 640 samples
  mouse <- mouse[, mouse@meta.data$class == "Glutamatergic"] # 7366 samples
}

mouse <- clustering_pipeline(mouse, clustering_resolution = 0.4) 
turtle <- clustering_pipeline(turtle, clustering_resolution=0.4) 
# Which groups/clusters to compare?
Idents(turtle) <- turtle@meta.data$seurat_clusters 
Idents(mouse) <- mouse@meta.data$subclass

# 1. For each dataset: marker genes findAllMarkers()
# To do: courser clustering of turtle neurons? Now have 17. 
# Or finer clustering of mouse neurons. 
top.mouse.markers <- get_top_markers(mouse) # Gaba: 2588, Glut: 4798
top.turtle.markers <- get_top_markers(turtle) # 2329, 2343

# 3. Intersect to select shared marker genes
shared.markers <- intersect(top.mouse.markers, top.turtle.markers)
print(length(shared.markers)) # ~419, 636

# 4. Average across clusters. slot = 'counts' -> don't exponentiate
# 5. Log(x+1) + 0.1
# 6. : across-cluster mean
turtle.avg <- score_genes(turtle, top_markers = shared.markers)
mouse.avg <- score_genes(mouse, top_markers = shared.markers)

# 6. Correlate. 
C <- cor(turtle.avg, mouse.avg) # 18 x 7
heatmap(C)

# Fixed color range
# https://stackoverflow.com/questions/17820143/how-to-change-heatmap-2-color-range-in-r
colors = seq(-.4,.4,length=100)
heatmap.2(as.matrix(C), col = hcl.colors(99, "RdYlBu",rev=TRUE), 
         breaks=colors, density.info="none", trace="none", 
         dendrogram=c("row"), symm=F,symkey=F,symbreaks=T)


# 7. Redo after shuffling to get significance. 
# To do: mouse - merge serpinf1 into Vip. 

