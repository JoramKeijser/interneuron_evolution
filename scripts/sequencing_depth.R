## Investigate sequencing depth
library(dplyr)
library(ggplot2)
library(Seurat)
library(SeuratDisk)
library(data.table)

datadir <- "./data/"
# Load data
mouse <- LoadH5Seurat(paste0(datadir, "Tasic18/", "mouse.h5seurat"))
mouse # 42094 genes x 15413 cells

# count depth aka n_c: total number of counts per cell
nc <- colSums(mouse)
plot(nc, mouse@meta.data$total_reads)
hist(nc)
median(nc) #1,635,968 = 1.6e6
# count fraction (p_gc): fraction of counts in each cell
p_g <- rowSums(mouse) / sum(nc)
hist(p_g)
median(p_g) # 2e-7

# expected value
mu_cg = outer(p_g, nc) # genes x cells
dim(mu_cg)
hist(mu_cg)

receptor_genes <- c("CHRM3", "CHRM4", "CHRNA4", "CHRNA5")
p_g[receptor_genes]

mu_chr <- outer(nc, p_g[receptor_genes])
mu_chr <- as.data.frame(mu_chr)
hist(mu_chr$CHRM4)

true_counts <- GetAssayData(object = mouse, slot = "counts")[receptor_genes, ]
#true_counts <- t(true_counts)
true_counts <- as.data.frame(true_counts)
true_counts <- transpose(true_counts)
colnames(true_counts) <- receptor_genes
hist(true_counts$CHRM3, breaks=1000)

# Now load the bird data, see how many counts we'd expct
bird <- LoadH5Seurat("./data/Colquitt21/bird.h5seurat") 
bird # 218342 features across 29497 samples within 1 assay 
# Select largest dataset
#bird <- bird[, bird@meta.data$species == 'zf']
nc_bird <- colSums(bird)
mu_chr_bird <- outer(nc_bird, p_g[receptor_genes])
mu_chr_bird <- as.data.frame(mu_chr_bird)
colnames(mu_chr_bird )
hist(mu_chr_bird$CHRM3)

# Turtle
turtle <- LoadH5Seurat("./data/Tosches18_turtle/turtle.h5seurat") 
turtle <- turtle[, turtle@meta.data$clusters %>% startsWith("i")]  # 640 samples
nc_turtle <- colSums(turtle)
mu_chr_turtle <- outer(nc_turtle, p_g[receptor_genes])
mu_chr_turtle <- as.data.frame(mu_chr_turtle)
hist(mu_chr_turtle$CHRM4)
true_counts_turtle <- GetAssayData(object = turtle, slot = "counts")[receptor_genes, ]
true_counts_turtle <- as.data.frame(true_counts_turtle)

# Human
source("./scripts/human2seurat.R")
nc_human <- colSums(human)
mu_chr_human <- outer(nc_human, p_g[receptor_genes])
mu_chr_human <- as.data.frame(mu_chr_human)
hist(mu_chr_human$CHRM3) # More counts

# Now generate poisson realizations with these means
generate_poission_counts <- function(x, gene_name, organism){
  counts <- rpois(n=dim(x)[1], lambda = x[[gene_name]])
  hist(counts, breaks=100, main = paste0(organism, ": ", gene_name))
  #print(table(counts)) # 24%   
  print(paste(organism, ": ", gene_name, round(mean(counts > 0)*100), "%"))
}

for (gene in receptor_genes){
  generate_poission_counts(mu_chr, gene, "mouse")
  generate_poission_counts(mu_chr_human, gene, "human")
  generate_poission_counts(mu_chr_bird, gene, "bird")
  generate_poission_counts(mu_chr_turtle, gene, "turtle")
  print("")
}

# But: overestimate in mouse. each receptor expressed by 99% of the neurons
# Need to be more specific, and do it for each cell type? 
# Or zero-inflated model?
# Shouldn't smear out counts over all neurons
round(colMeans(true_counts>0) * 100)
#  CHRM3  CHRM4 CHRNA4 CHRNA5 
#   70     25     61      9 
hist(true_counts$CHRM3)
hist(true_counts[true_counts$CHRM3>0,]$CHRM3)

true_counts_human <- GetAssayData(object = human, slot = "counts")[receptor_genes, ]
round(rowMeans(true_counts_human>0) * 100)
hist(true_counts_human['CHRM3',])

# CHRM3  CHRM4 CHRNA4 CHRNA5 
# 72      0     13      2 
# CHRM4 really does seem absent in humans

true_counts_turtle <- GetAssayData(object = turtle, slot = "counts")[receptor_genes, ]
round(rowMeans(true_counts_turtle>0) * 100)
# CHRM4 is actually more abundant 
hist(true_counts_turtle['CHRM3',])


true_counts_bird <- GetAssayData(object = bird, slot = "counts")[receptor_genes, ]
round(rowMeans(true_counts_bird>0) * 100)
# CHRM4 is actually more abundant 
hist(true_counts_bird['CHRM3',])

# (1) Extend model to cell type?
# (2) Hypothesis testing. Is Chrm4 count really more than expected from mouse?

