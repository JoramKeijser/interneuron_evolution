#!/bin/bash
# Assumes we're in the repo directory 

# Create subdirectory for data
mkdir data
mkdir data/raw data/seurat data/anndata
# Next, download raw datasets into /raw
cd data/raw

# Mouse data from Tasic et al. 
echo "Download dataset 1/4: mouse"
mkdir tasic
cd tasic
wget https://celltypes.brain-map.org/api/v2/well_known_file_download/694413985
file-roller -h 694413985 # extract
mv 694413985_FILES/* ./ # move subfolder -> folder
rm -r 694413985_FILES 694413985	 # clean up

# Human data from Bakken et al. 
echo "Download dataset 2/4: human"
mkdir ../bakken
cd ../bakken
wget https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/metadata.csv
wget https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/matrix.csv

# Bird data from Colquitt et al. 
echo "Download dataset 3/4: songbird"
mkdir ../colquitt
cd ../colquit
wget -O HVC_RA_RNA_counts.csv 'https://cloud.biohpc.swmed.edu/index.php/s/nLicEtkmjGGmRF8/download?path=%2FHVC_RA&files=HVC_RA_RNA_counts.csv&downloadStartSecret=6cncz7uwc0g'

# Turtle data from Tosches et al. 
echo "Download dataset 4/4: turtle"
mkdir ../tosches
cd ../tosches
wget https://public.brain.mpg.de/Laurent/ReptilePallium2018/turtle.neurons.Robj

cd ../../.. # return to main dir. 

