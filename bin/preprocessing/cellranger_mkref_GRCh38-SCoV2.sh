#!/bin/bash

# set parameters for SGE submission
#$ -N mkfastq
#$ -l arch=linux-x64
#$ -pe multislot 20
#$ -b n
#$ -q all.q
#$ -o log/mkref.log
#$ -e log/mkref.log
#$ -cwd

# Choose PATH for the cellranger version
export PATH=/vol/biotools/bin:$PATH
export PATH=/vol/projects/odietric/bin/cellranger/cellranger-5.0.1:$PATH

# Show at beginning of job
echo ""
hostname -f
echo "Starting on" $(date)
echo ""

# Set variables
fa_1=refdata-gex-GRCh38-2020-A/fasta/genome.fa
gtf_1=refdata-gex-GRCh38-2020-A/genes/genes.gtf

fa_2=SCoV2.fa
gtf_2=SCoV2.gtf

fa=GRCh38-SCoV2.fa
gtf=GRCh38-SCoV2.gtf
out=GRCh38-SCoV2

# Go to output directory
cd data/genomes

# Combine fasta and gtf files
cat $fa_1 > $fa
cat $fa_2 >> $fa
# grep ">" $fa checks for FASTA entries

cat $gtf_1 > $gtf
cat $gtf_2 >> $gtf

# Create reference
cellranger mkref --genome=$out --fasta=$fa --genes=$gtf --memgb=32

echo ""
echo "Done on" $(date)
