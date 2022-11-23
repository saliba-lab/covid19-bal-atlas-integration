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
raw=data/BCB/raw.h5ad
filtered=data/BCB/cohortA.h5ad

# Setup
python bin/method_filter.py -o $filtered --cohortA $raw
python bin/method_normalization.py $filtered
python bin/method_hvfeatures.py $filtered

# Integration
python bin/method_integration-PCA.py $filtered
Rscript bin/method_integration-fastMNN.R $filtered
python bin/method_integration-scVI.py $filtered

# Mapping
export PATH=~/miniconda3/envs/covid19-bal-atlas-scarches/bin:$PATH
bash bin/dataset-HLCA-core.sh
python bin/method_mapping-HLCA.py $filtered

# Annotation
export PATH=~/miniconda3/envs/covid19-bal-atlas/bin:$PATH
Rscript bin/method_embed-cluster.R -o $plot_ann $core

# Metrics
export PATH=~/miniconda3/envs/covid19-bal-atlas-scib/bin:$PATH
python bin/metric_scIB.py $filtered
