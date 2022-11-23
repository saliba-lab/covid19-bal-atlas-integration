#!/bin/bash

#SBATCH --job-name=count-viral
#SBATCH --output=log/cellranger/count-viral.log
#SBATCH --error=log/cellranger/count-viral.log
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=48:00:00		# walltime (HH:MM:SS)
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oliver.dietrich@helmholtz-hiri.de
#SBATCH --clusters=bioinf

# Choose PATH for the cellranger version
export PATH=/vol/biotools/bin:$PATH
export PATH=/vol/projects/odietric/bin/cellranger/cellranger-5.0.1:$PATH

# Set variables
base_dir=$(pwd)
fastq_path=$base_dir/docs/fastq-path.csv
out=$base_dir/data/processed/viral/
ref=$base_dir/data/genomes/viral-assembly

# Create output dir
mkdir $out

# Move to output directory
cd $out

# Run cellranger for each run
for i in $(cat $fastq_path)
do

  # User message
  echo $(date)

  # Checkpoint --- Samplesheets and FASTQs
  i=$base_dir/$i
  if [[ -d $i ]]; then
    sample=$(basename $i)
    fastqs=$i
    echo "Sample $sample and fastqs ($fastqs) are present"
    trigger=1
  fi

  # Checkpoint --- Count matrix ---
  if [[ -e $sample/outs/molecule_info.h5 ]]; then
    echo "Output already present. Skipping sample $sample"
    trigger=0
  elif [[ -d $sample ]]; then
    echo "Output directory present but pipeline did not finish successfully. Removing directory."
    rm -r $sample
  fi

  # Run cellranger
  if [[ $trigger == 1 ]]; then
    cellranger count --id=$sample --fastqs=$fastqs --sample=$sample --transcriptome=$ref --nosecondary
  fi

done
