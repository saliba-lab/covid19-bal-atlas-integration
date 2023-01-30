#!/bin/bash

#SBATCH --job-name=jupyter
#SBATCH --output=log/jupyter.log
#SBATCH --error=log/jupyter.log
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=10:00:00         # walltime (HH:MM:SS)
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oliver.dietrich@helmholtz-hiri.de
#SBATCH --clusters=bioinf

# Set up variables
env=covid19-bal-atlas-integration
conda_path=/home/$USER/miniconda3

# Activate local environment
unset PYTHONPATH
export PATH=$conda_path/envs/$env/bin:$PATH

# Run jupyter notebook
jupyter notebook --no-browser --port=8080
