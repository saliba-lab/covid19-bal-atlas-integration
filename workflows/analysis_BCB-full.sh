#!/bin/bash

#SBATCH --job-name=analysis_BCB-full
#SBATCH --output=log/analysis_BCB-full.log
#SBATCH --error=log/analysis_BCB-full.log
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

unset PYTHONPATH

# Variables
raw=data/BCB/raw.h5ad
csv=analysis/BCB/qc/colData.csv
filtered=data/BCB/full.h5ad

unset PYTHONPATH

# Setup
export PATH=~/miniconda3/envs/covid19-bal-atlas-integration/bin:$PATH
python bin/method_filter.py -o $filtered $raw
python bin/method_normalization.py $filtered
python bin/method_hvfeatures.py $filtered

# Integration
python bin/method_integration-PCA.py $filtered
python bin/method_integration-scVI.py $filtered # use --gpu if present

# Integration (R)
export PATH=~/miniconda3/envs/covid19-bal-atlas-scran/bin:$PATH
Rscript bin/method_integration-fastMNN.R $filtered

# Mapping
export PATH=~/miniconda3/envs/covid19-bal-atlas-scarches/bin:$PATH
bash bin/dataset-HLCA-core.sh
python bin/method_mapping-HLCA.py $filtered

# Annotation
export PATH=~/miniconda3/envs/covid19-bal-atlas-scran/bin:$PATH
Rscript bin/method_scDblFinder.R $filtered
Rscript bin/method_embed-cluster.R -o $plot_ann $core

# Metrics
export PATH=~/miniconda3/envs/covid19-bal-atlas-scib/bin:$PATH
python bin/metric_scIB.py $filtered

