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
1. Install conda (follow instructions and accept defaults)
      ```
      curl -o miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
      bash miniconda.sh
      rm miniconda.sh
      ```
3. Create conda environments
   ```
   conda env create -f envs/default.yml
   ```
