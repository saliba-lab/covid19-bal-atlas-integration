#!/bin/bash

#SBATCH --job-name=count
#SBATCH --output=log/cellranger_count.log
#SBATCH --error=log/cellranger_count.log
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=48:00:00         # walltime (HH:MM:SS)
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oliver.dietrich@helmholtz-hiri.de
#SBATCH --clusters=bioinf

# Choose PATH for the cellranger version
export PATH=/vol/biotools/bin:$PATH
export PATH=/vol/biotools/bin/cellranger-7.0.0:$PATH

# Define variables
base_dir=$(pwd)
samplesheets=$base_dir/data/samplesheets/count
out=$base_dir/data/libraries/
ref=$base_dir/data/genomes/GRCh38-viral
htoref=/home/odietric/SIGA/databases/cellranger_refs/HTO-features.csv

# Create output directory
if [[ -e $out ]]; then
  echo "Output directory: $out"
else
  echo "Creating output directory in $out"
  mkdir $out
fi

# Make sure samplesheets are up to date
echo "Create samplesheets from overview"
Rscript reports/overview_BCB.R

# Move to output directory
cd $out

# Run cellranger for each sample
for i in $samplesheets/*
do
  echo $(date)

  sample=$(basename $i)
  sample=${sample%.csv}
  echo $sample
  trigger=1

  # Checkpoint
  # Presence of output files
  if [[ -e $sample/outs/molecule_info.h5 ]]; then
    echo "Output already present. Skipping sample $sample"
    trigger=0
  elif [[ -d $sample ]]; then
    echo "Output directory present but pipeline did not finish successfully. Removing directory."
    rm -r $sample
  fi

  if [[ $trigger == 1 ]]; then
    cellranger count --id=$sample --libraries=$i --transcriptome=$ref --feature-ref=$htoref \
                     --nosecondary
  fi


done
