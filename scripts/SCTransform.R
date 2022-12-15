library(dplyr)
library(Seurat)
library(SeuratDisk)

# Independently for each sample? 
# https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1

bird <- LoadH5Seurat("./data/raw/bird/bird.h5seurat")
# SCTransform 
bird <- SCTransform(bird, vars.to.regress = "percent.mito", 
                    verbose = TRUE, variable.features.n = 3000)
SaveH5Seurat(bird, "./data/raw/bird/birdSCT")

# mouse
mouse <- LoadH5Seurat("./data/raw/mouse/mouse.h5seurat")
colnames(mouse@meta.data)
mouse <- SCTransform(mouse, vars.to.regress = "percent_mt_exon_reads", 
                     verbose=TRUE, variable.features.n = 3000)
mouse
SaveH5Seurat(mouse, "./data/raw/mouse/mouseSCT")

# To do: turtle. 
