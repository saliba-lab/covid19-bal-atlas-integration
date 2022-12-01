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

# Choose PATH for the cellranger version
export PATH=/vol/biotools/bin:$PATH
export PATH=/vol/projects/odietric/bin/cellranger/cellranger-5.0.1:$PATH

# Set variables
base_dir=$(pwd)
out=$base_dir/data/processed
csv=$base_dir/docs/combined_aggr.csv
name=combined

# Create aggr.csv from all samples
echo "library_id,molecule_h5" > $csv
for i in $out/*
do

  sample=$(basename $i)
  file=$i/outs/molecule_info.h5
  if [[ -e $file ]]; then
    echo $sample,$file >> $csv
  fi

done

# Move to output directory
cd $out

# Checkpoint
if [[ -d $name ]]; then
  echo "Output directory already exists. Pipeline will not be executed."
  trigger=0
else
  trigger=1
fi

# Run cellranger aggr
if [[ $trigger == 1 ]]; then
  cellranger aggr --id=$name --csv=$csv --normalize=none --nosecondary
fi
