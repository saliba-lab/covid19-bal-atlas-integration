# Atlas-level data integration of BAL samples from severe COVID-19

This repository contains code for the processing and analysis of scRNA-seq samples from patients with severe COVID-19.

## Directory Structure

* `bin/` - Scripts used for the analysis
* `docs/` - Document files
* `data/` - Data output (large, not included in **git**)
* `analysis/` - Analysis output (not included in **git**)
* `envs/` - **conda** environment YAML files
* `LICENSE` - The project license

## Analysis workflow

1. Create the conda environment
   ```
   conda env create -f envs/_template.yml
   ```
