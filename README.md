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
   python bin/dataset-combined.py -f data/combined.h5ad
   Rscript bin/analysis-combined-qc.R
   python bin/analysis-combined-qc.py
   python bin/analysis-combined-scVI.py
   Rscript bin/analysis-combined-fastMNN.R
   Rscript bin/analysis-combined-celltype.R
   python bin/analysis-combined-explore.py
   ```
1. Mapping to the Human Lung Cell Atlas
   ```
   bash bin/dataset-hlca.sh
   ```
## Reference and Query Datasets

Different datasets will be used for integration and mapping. Information about the datasets can be found in the [issues](https://github.com/saliba-lab/covid19-bal-atlas-integration/issues) while the download and use of the datasets is performed using scripts in `bin/` such as
> [dataset-combined.py](https://github.com/OliverDietrich/covid19.atlas/blob/main/bin/dataset-combined.py)

## Raw data processing

Raw data in this case means Illumina base call (bcl) files. The processing is performed by cellranger pipelines run on a [Slurm](https://slurm.schedmd.com/overview.html) controlled high performance cluster (HPC).

### Combined BAL dataset

1. Create genome reference
   ```
   sbatch bin/GRCh38-SCoV2_cellranger-mkref.sh
   ```
1. Run cellranger pipelines
   ```
   sbatch bin/combined_cellranger-mkfastq.sh
   sbatch bin/combined_cellranger-count.sh
   sbatch bin/combined_cellranger-aggr.sh
   ```
1. Share data
   Take data matrix in data/processed/combined/outs/filtered_feature_bc_matrix.h5 and upload it to cloud service (we use [Nubes](https://nubes.helmholtz-berlin.de))
   
### Viral assembly

1. Create genome reference
   ```
   Rscript bin/download-viral-genomes.R
   sbatch bin/GRCh38_SCoV2_cellranger-mkref.sh
   ```
1. Run cellranger pipelines
   ```
   sbatch bin/combined-viral_cellranger-mkfastq.sh
   sbatch bin/combined-viral_cellranger-count.sh
   sbatch bin/combined-viral_cellranger-aggr.sh
   ```
1. Share data
   Take data matrix in data/processed/combined/outs/filtered_feature_bc_matrix.h5 and upload it to cloud service (we use [Nubes](https://nubes.helmholtz-berlin.de))
