#!/bin/bash
# Preprocess all the datasets: create Seurat and AnnData files
# Assumes we're in the repo directory and raw data is in data/raw

echo "Mouse"
Rscript scripts/preprocess_mouse.R

echo "Human"
Rscript scripts/preprocess_human.R

echo "Bird"
Rscript scripts/preprocess_bird.R

echo "Turtle"
Rscript scripts/preprocess_turtle.R

