# Atlas-level data integration of BAL samples from severe COVID-19

This repository contains code for the processing and analysis of scRNA-seq samples from patients with severe COVID-19.

## Directory Structure

* `bin/` - Scripts used for the analysis
* `docs/` - Document files
* `data/` - Data output (large, not included in **git**)
* `workflows/` - Workflows executing scripts in **bin/**
* `reports/` - Reports for exploratory analysis
* `analysis/` - Analysis output (not included in **git**)
* `man/` - Manuals explaining the rationale behind analyses
* `envs/` - **conda** environment YAML files
* `LICENSE` - The project license

## Data Availability

Shared on web portal: TODO

Download as part of workflows, described in bin/dataset_*

## Analysis workflows

The analysis workflow include download & setup of datasets. Quality control, data integration, cell annotation, 
evaluation using (scIB) metrics and visualization of reports.

The workflows are written in bash and can be submitted using SLURM.

To reproduce the analysis please
1. Install [conda](https://docs.conda.io/en/latest/miniconda.html#) (follow instructions and accept defaults), close and re-open the console
   ```
   curl -o miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash miniconda.sh
   rm miniconda.sh
   ```
1. Clone and enter the git repository
   ```
   git clone https://github.com/saliba-lab/covid19-bal-atlas-integration.git
   cd covid19-bal-atlas-integration
   ```
1. Create conda environments
   ```
   conda env create -f envs/default.yml
   conda activate covid19-bal-atlas
   mamba init
   ```
   close and re-open terminal
   ```
   mamba env create -f envs/scran.yml
   mamba env create -f envs/scArches.yml
   mamba env create -f envs/scIB.yml
   ```
   
## Addendum

To use nextflow workflows please install java
   ```
   sudo apt install default-jre
   ```
