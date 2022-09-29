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

Different datasets will be used for integration and mapping. Information about the datasets can be found here while the download and use of the datasets is performed using scripts in `bin/` such as
> dataset-combined.py

### References

Reference datases are well annotated atlases of the healthy or diseased lung

#### HLCA

The Integrated Human Lung Cell Atlas (HLCA). Read up on the [GitHub](https://github.com/LungCellAtlas/HLCA) and [preprint](https://www.biorxiv.org/content/10.1101/2022.03.10.483747v1).

1. HLCA core:
   > https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293
1. HLCA extended
   > https://beta.fastgenomics.org/datasets/detail-dataset-427f1eee6dd44f50bae1ab13f0f3c6a9#Files

### Querys

## Raw data processing

Raw data in this case means Illumina base call (bcl) files. The processing is performed by cellranger pipelines run on a [Slurm](https://slurm.schedmd.com/overview.html) controlled high performance cluster (HPC).
