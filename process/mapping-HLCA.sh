#!/bin/bash

#SBATCH --job-name=mapping-HLCA
#SBATCH --output=log/mapping-HLCA.log
#SBATCH --error=log/mapping-HLCA.log
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=48:00:00         # walltime (HH:MM:SS)
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oliver.dietrich@helmholtz-hiri.de
#SBATCH --clusters=bioinf

# Variables
file=$1

# Set PATHs
unset PYTHONPATH
export PATH=$HOME/miniconda3/envs/covid19-bal-atlas-scarches/bin:$PATH

# Download reference data
bash bin/dataset-HLCA-core.sh

# Run mapping
python bin/method_mapping-HLCA.py $file
