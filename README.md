## Interneuron evolution
Analysing the evolution of cortical interneurons using single cell RNA sequencing data. This code accompanies a forthcoming paper in which we constrast an functional (or optimisation) view of interneurons with evolutionary-developmental data. 

![figure1](./figures/paper_figs/figure1.png)

## Installation

Make a copy of this repo cd into the root folder of the repo, and download the raw data.
```
git clone https://github.com/JoramKeijser/elfn1_evolution/
cd elfn1_evolution
./scripts/download_data.sh
```
Recreate the conda environment to install the required Python and R packages/libraries. 
```
conda env create --name elfn1_env --file environment.yml
```
Install the project:
```
pip install -e .
```

## Organization
The code is organized into the following folders:

- `scripts` contains scripts that run the analyses
- `src` contains the code that is shared by several scripts
- `data` will contain the data after running download_data.sh (see Installation above)
- `figures` contains the figures
- `results` contains text files 
 
## Analysis 

The analysis pipeline is shown in the figure below. Each step corresponds to one or several scripts. 

![figures1](./figures/paper_figs/figures1.png)

Execute the following code to run the pipeline from beginning to end:
```
cd elfn1_evolution
conda activate elfn1_env
./scripts/download_data.sh # Downloads data into data/raw
./scrtipts/preprocess.sh # Puts data into Seurat and AnnData files
python ./scripts/visualize_expression.py # UMAPs and violin plots
Rscript ./scripts/correlation.R # correlate gene expression of birds and mice
Rscript ./scripts/integrate.R # integrate bird and mouse data
python ./scripts/integrated_umaps.py # visualize integrated data
```

## Acknowledgements

Thanks to Patrick Mineault for writing the [Good Research Codebook](https://goodresearch.dev/) that we used as guideline for developing this code. 
