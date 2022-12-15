# Save raw bird data as Seurat
# https://satijalab.org/seurat/articles/integration_introduction.html
library(dplyr)
library(Seurat)
library(SeuratDisk)

rm(list = ls())
load("./data/raw/turtle/turtle.neurons.Robj", verbose=TRUE)
    
data.attributes <- attributes(turtle.neurons)
# An old seurat object
# 21167 genes across 5901 samples
counts <- attr(turtle.neurons, "raw.data") # genes x cells
cell.types <- attr(turtle.neurons, "ident")
cell.names <- attr(turtle.neurons, "cell.names")
gene.names <- rownames("counts") #obj.attributes[["counts"]]@Dimnames[[1]]
data.info <- attr(turtle.neurons, "data.info")

data.dir <- "./data/raw/turtle/"
write.table(data.info, "./data/raw/turtle/metadata.csv", sep=",")
write.table(counts, "./data/raw/turtle/counts.csv", sep=",")

turtle <- CreateSeuratObject(counts = counts, project = "turtle", 
                             meta.data = data.info, 
                             min.cells = 1, min.features = 1)
turtle
SaveH5Seurat(turtle, "./data/raw/turtle/turtle")


