#!/bin/bash

#SBATCH --job-name=setup
#SBATCH --output=log/workflow/setup.log
#SBATCH --error=log/workflow/setup.log
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
csv=data/BCh/qc_colData.csv
plot_qc=analysis/BCB/qc

# Download data
export PATH=~/miniconda3/envs/covid19-bal-atlas-scib/bin:$PATH
python bin/dataset-BCB.py $raw

# Calculate QC metrics
export PATH=~/miniconda3/envs/covid19-bal-atlas-scran/bin:$PATH
Rscript bin/method_qc.R $raw
