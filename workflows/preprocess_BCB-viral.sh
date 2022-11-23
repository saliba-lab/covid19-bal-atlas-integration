#!/bin/bash

#SBATCH --job-name=preprocess_BCB-viral
#SBATCH --output=log/preprocess/BCB-viral.log
#SBATCH --error=log/preprocess/BCB-viral.log
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

Rscript bin/download-viral-genomes.R
bash bin/GRCh38_SCoV2_cellranger-mkref.sh
bash bin/combined-viral_cellranger-mkfastq.sh
bash bin/combined-viral_cellranger-count.sh
bash bin/combined-viral_cellranger-aggr.sh
