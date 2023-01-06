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
Recreate the conda environment:
```
conda env create -n cb --file environment.yml
```
Install the project:
```
pip install -e .
```

## Organization
The code is organized into the following folders:

- `scripts` contains scripts that run the analyses
- `src` contains the code that is shared by several sripts
- `data` will contain the data after running the code under `Installation`
- `figures` contains the figures
- `results` contains text files 
 
## Analysis 

The analysis pipeline is shown in the figure below. Each step corresponds to one or several scripts. 

![figures1](./figures/paper_figs/figures1.png)

Preprocess the raw data by putting them into AnnData and Seurat files. Then run the other functions
