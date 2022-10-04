#!/bin/bash

#SBATCH --job-name=mkref-viral
#SBATCH --output=log/cellranger/mkref-viral.log
#SBATCH --error=log/cellranger/mkref-viral.log
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=24:00:00		# walltime (HH:MM:SS)
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oliver.dietrich@helmholtz-hiri.de
#SBATCH --clusters=bioinf

# Choose PATH for the cellranger version
export PATH=/vol/biotools/bin:$PATH
export PATH=/vol/projects/odietric/bin/cellranger/cellranger-5.0.1:$PATH

# Show at beginning of job
hostname
base_dir=$(pwd)

# Convert gff3 to gtf
for i in $base_dir/data/genomes/viral/*.gff3
do
gffread $i -T -o ${i/gff3/gtf}
done

# Go to output directory
cd $base_dir/data/genomes

# Output files
out=viral-assembly

# Input files
genomes=$base_dir/data/genomes/viral
fasta=$genomes/viral-assembly.fasta
rm $fasta
gtf=$genomes/viral-assembly.gtf # manually curated from converted gtf files
cat $genomes/*.fasta > $fasta # pasted from individual fasta files

# Create reference
cellranger mkref --genome=$out --fasta=$fasta --genes=$gtf --memgb=32
