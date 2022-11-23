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

# Set PATH
export PATH=~/miniconda3/envs/covid19-bal-atlas/bin:$PATH

# Variables
raw=data/BCh/raw.h5ad
csv=data/BCh/qc_colData.csv
filtered=data/BCh/full.h5ad

# Run scripts
python bin/method_filter.py -o $filtered --csv $csv $raw
