## Interneuron evolution
Analysing the evolution of cortical interneurons using single cell RNA sequencing data. This code accompanies a forthcoming paper in which we constrast an functional (or optimisation) view of interneurons with evolutionary-developmental data. 

![figure1](./figures/paper_figs/figure1.png)

## Installation

Clone this repository:
```
git clone https://github.com/JoramKeijser/elfn1_evolution/
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

Cd to the repo directory, and download the data:
```
cd elfn1_evolution
./scripts/download_data.sh
```
Run the entire pipeline:
```
./run_pipeline.sh
``` 
