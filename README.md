# Atlas-level data integration of BAL samples from severe COVID-19

This repository contains code for the processing and analysis of scRNA-seq samples from patients with severe COVID-19.

## Directory Structure

* `bin/` - Scripts used for the analysis
* `docs/` - Document files
* `data/` - Data output (large, not included in **git**)
* `analysis/` - Analysis output (not included in **git**)
* `man/` - Manuals explaining the rationale behind analyses
* `envs/` - **conda** environment YAML files
* `LICENSE` - The project license

## Analysis workflow

The analyis workflow starts from the cell x gene count matrix.

1. Create the conda environment
   ```
   conda env create -f envs/_template.yml
   conda activate covid19-bal-atlas
   ```
1. Run the scripts for analysis of the combined BAL dataset
   ```
   python bin/dataset-combined.py
   ```

## Reference and Query Datasets

Different datasets will be used for integration and mapping. Information about the datasets can be found in the [issues](https://github.com/saliba-lab/covid19-bal-atlas-integration/issues) while the download and use of the datasets is performed using scripts in `bin/` such as
> [dataset-combined.py](https://github.com/OliverDietrich/covid19.atlas/blob/main/bin/dataset-combined.py)

## Raw data processing

Raw data in this case means Illumina base call (bcl) files. The processing is performed by cellranger pipelines run on a [Slurm](https://slurm.schedmd.com/overview.html) controlled high performance cluster (HPC).
