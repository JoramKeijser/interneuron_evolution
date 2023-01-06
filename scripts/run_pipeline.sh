#!/bin/bash
# Run the entire pipeline

# First download the data
cd elfn1_evolution # assumes you're in the base dir
./download_data.sh

# Run the actual code
conda activate elfn1_env 
./scripts/preprocess.sh # Puts data into Seurat and AnnData files
python ./scripts/visualize_expression.py # UMAPs and violin plots
Rscript ./scripts/correlate.R # correlate gene expression of birds and mice
Rscript ./scripts/integrate.R # integrate bird and mouse data
python ./scripts/visualize_integration.py # visualize integrated data

