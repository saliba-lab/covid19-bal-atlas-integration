#!/bin/bash

#SBATCH --job-name=analysis_combined
#SBATCH --output=log/analysis_combined.log
#SBATCH --error=log/analysis_combined.logs
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

# Variables
raw=data/BCB/raw.h5ad
csv=data/BCB/qc_colData.csv
filtered=data/BCB/full.h5ad

# Setup
export PATH=~/miniconda3/envs/covid19-bal-atlas-integration/bin:$PATH
python bin/method_filter.py -o $filtered $raw
python bin/method_normalization.py $filtered
python bin/method_hvfeatures.py $filtered

# Integration
python bin/method_integration-PCA.py $filtered
python bin/method_integration-scVI.py --gpu $filtered # use --gpu if present

# Integration (R)
export PATH=~/miniconda3/envs/covid19-bal-atlas-scran/bin:$PATH
Rscript bin/method_integration-fastMNN.R $filtered

# Mapping
export PATH=~/miniconda3/envs/covid19-bal-atlas-scarches/bin:$PATH
bash bin/dataset-HLCA-core.sh
python bin/method_mapping-HLCA.py $filtered

# Annotation
export PATH=~/miniconda3/envs/covid19-bal-atlas-scran/bin:$PATH
Rscript bin/method_doubletFinder.R $filtered
Rscript bin/method_embed-cluster.R -o $plot_ann $filtered

# Metrics
export PATH=~/miniconda3/envs/covid19-bal-atlas-scib/bin:$PATH
python bin/metric_scIB.py $filtered

