#!/bin/bash

#SBATCH --job-name=aggr
#SBATCH --output=log/cellranger/aggr.log
#SBATCH --error=log/cellranger/aggr.log
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=48:00:00         # walltime (HH:MM:SS)
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oliver.dietrich@helmholtz-hiri.de
#SBATCH --clusters=bioinf

# Set PATH
export PATH=~/miniconda3/envs/covid19-bal-atlas/bin:$PATH

# Variables
raw=data/BCh/raw.h5ad
csv=data/BCh/qc_colData.csv
filtered=data/BCh/initial.h5ad
plot_ann=analysis/BCh/initial/celltypes

# Setup
python bin/method_filter.py -o $filtered --csv $csv --subset $raw
python bin/method_normalization.py $filtered
python bin/method_hvfeatures.py $filtered

# Integration
python bin/method_integration-PCA.py $filtered
Rscript bin/method_integration-fastMNN.R $filtered
python bin/method_integration-scVI.py $filtered

# Mapping
export PATH=~/miniconda3/envs/covid19-bal-atlas-scarches/bin:$PATH
python bin/method_mapping-HLCA.py $filtered

# Metrics
export PATH=~/miniconda3/envs/covid19-bal-atlas-scib/bin:$PATH
python bin/method_integration-metrics.py $mapped

# Annotation
# Rscript bin/method_embed-cluster.R -o $plot_ann $core
