#!/bin/bash

#SBATCH --job-name=mkref
#SBATCH --output=log/cellranger_mkref.log
#SBATCH --error=log/cellranger_mkref.log
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=48:00:00         # walltime (HH:MM:SS)
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oliver.dietrich@helmholtz-hiri.de
#SBATCH --clusters=bioinf

# Choose PATH for the cellranger version
export PATH=/vol/biotools/bin:$PATH
export PATH=/vol/biotools/bin/cellranger-7.0.0:$PATH

# Show at beginning of job
hostname
echo "Starting on" $(date)

# Set variables
base_dir=$(pwd)
ref_1=/home/odietric/SIGA/databases/cellranger_refs/refdata-gex-GRCh38-2020-A/
fa_1=$ref_1/fasta/genome.fa
gtf_1=$ref_1/genes/genes.gtf

fa_2=$base_dir/docs/genomes/viral.fasta
gtf_2=$base_dir/docs/genomes/viral.gtf

fa=GRCh38-viral.fa
gtf=GRCh38-viral.gtf
out=GRCh38-viral

# Go to output directory
cd data/genomes

# FASTA
if [[ -e $fa ]]; then
  echo "FASTA present."
else
  echo "Combining FASTA files..."
  cat $fa_1 > $fa
  cat $fa_2 >> $fa
fi
# grep ">" $fa checks for FASTA entries

# GTF
if [[ -e $gtf ]]; then
  echo "GTF present."
else
  echo "Combining GTF files..."
  cat $gtf_1 > $gtf
  cat $gtf_2 >> $gtf
fi

# Make reference genome
if [[ -d $out ]]; then
  echo "$out already exists. Exiting."
  exit 1
else
  echo "Creating reference..."
  cellranger mkref --genome=$out --fasta=$fa --genes=$gtf --memgb=32
fi

echo "Done on" $(date)
