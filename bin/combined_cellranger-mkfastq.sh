#!/bin/bash

#SBATCH --job-name=mkfastq
#SBATCH --output=log/cellranger/mkfastq.log
#SBATCH --error=log/cellranger/mkfastq.log
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=48:00:00		# walltime (HH:MM:SS)
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oliver.dietrich@helmholtz-hiri.de
#SBATCH --clusters=bioinf

# Choose PATH for the cellranger version
export PATH=/vol/biotools/bin:$PATH
export PATH=/vol/projects/odietric/bin/cellranger/cellranger-5.0.1:$PATH

# Show at beginning of job
hostname -f
date

# Define file PATHs
basedir=$(pwd)
samplesheets="$basedir/docs/samplesheets/*"
bcl="$basedir/data/raw/bcl/"
fastq="data/raw/fastq/"

# Change directory (cellranger does not accept PATHs as --id argument)
cd $fastq
echo "Location for FASTQs:" $(pwd)

# Iterate for all runs that have samplesheets
for csv in $samplesheets
do
  # Set variables
  csv=$pwd$csv
  run=$(basename $csv)
  id=${run%.csv}
  run=$bcl$id

  echo ""
  echo "For id:" $id
  echo "samplesheet:" $csv

  # Checkpoint --- BCL files ---
  if [[ -d $run ]]; then
    echo $run "is present"
    trigger=1
  else
    echo "Warning:" $run "is not present"
    trigger=0
  fi

  # Checkpoint --- FASTQ files ---
  if [[ -d $id/outs/fastq_path ]]; then
    echo "FASTQ files are already present. $id will be skipped"
    trigger=0
  elif [[ -d $id ]]; then
    echo "FASTQ dir exists but pipeline did not finish successfully. Removing $id."
    rm -r $id
  fi

  # Checkpoint 3 --- Trigger ---
  if [[ $trigger == 1 ]]; then
    cellranger mkfastq --id=$id --run=$run --csv=$csv --delete-undetermined
  fi

done
