# Helper functions for analysis with Seurat
library(dplyr)
library(Seurat)
library(SeuratDisk)

# Standard workflow for clustering
clustering_pipeline <- function(seurat_obj, 
                                clustering_resolution = 0.8, umap_dims = 1:15){
  seurat_obj %>% NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>% 
    FindNeighbors() %>%
    FindClusters(resolution = clustering_resolution) %>%
    RunUMAP(dims=umap_dims)
}


get_cluster_top_markers <- function(cluster_name, marker_df){
  # Select top markers for particular cluster
  marker_df %>% filter(cluster == cluster_name) %>%
    filter(p_val_adj < 0.05) %>%
    arrange(desc(avg_log2FC)) %>% 
    head(400) %>%
    select(gene) 
}

get_top_markers <- function(seurat_obj){
  # Select top markers for all clusters
  marker_df <- FindAllMarkers(seurat_obj, min.pct=0.2, test.use='t', 
                              max.cells.per.ident = 200)
  cluster_names <- unique(Idents(seurat_obj))
  top.markers <- lapply(cluster_names, function(x){get_cluster_top_markers(x, marker_df)})
  unlist(unlist(top.markers, recursive = F), recursive=F)
}

# Specificity score
# 4. Average across clusters. slot = 'counts' -> don't exponentiate

score_genes <- function(seurat_obj, top_markers){
  # Average within Idents
  avg.expr <- AverageExpression(seurat_obj, features = top_markers,
                                slot='counts')$RNA
  avg.expr <- as.data.frame(avg.expr)
  # Log transform
  avg.expr <- log1p(avg.expr) + 0.1
  # Across-cluster avg
  avg.expr <- avg.expr / rowMeans(avg.expr)
  return(avg.expr)
}
