#!/bin/bash

#SBATCH --job-name=tar2bcl
#SBATCH --output=log/tar2bcl_221014_A00643_0567_AHYMN3DSX3_P1713.log
#SBATCH --error=log/tar2bcl_221014_A00643_0567_AHYMN3DSX3_P1713.log
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=48:00:00         # walltime (HH:MM:SS)
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oliver.dietrich@helmholtz-hiri.de
#SBATCH --clusters=bioinf

# Print hostname & date
hostname -f
date

# Global variables
in_dir="data/raw/tar/221014_A00643_0567_AHYMN3DSX3_P1713"
files=$(ls $in_dir/BCLS_LANE_*.tar)
out_dir="data/raw/bcl"
counter=""

if [[ -d $out_dir ]]; then
  echo "Extract files to $out_dir"
else
  mkdir $out_dir
fi

# Compare md5sums
for file in $files; do

  # Variables
  md5_true=$(cat $file.md5)
  md5_true=$(cut -d " " -f1 <<< $md5_true)
  md5_local=$(cat $file.md5loc)
  md5_local=$(cut -d " " -f1 <<<  $md5_local)

  echo "Unpacking $file ..."
  echo $md5_true
  echo $md5_local

  if [[ $md5_true == $md5_local ]]; then
    echo "md5sums check out. Adding '1' to counter."
    counter="$counter"1
  else
    echo "md5sums are different. Adding '0' to counter. Please download again."
    counter="$counter"0
  fi

done

if [[ -d $out_dir ]]; then
  echo "Directory $out_dir already exists. Exiting."
  exit 2
fi

if [[ $counter == 1111 ]]; then
  for file in $files; do
    tar -xvf $file --directory $out_dir
  done
else
  echo "Counter is $counter. Expecting 1111. Some samples need to be downloaded again."
fi
