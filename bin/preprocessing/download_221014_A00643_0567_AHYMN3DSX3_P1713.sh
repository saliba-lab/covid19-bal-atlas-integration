#!/bin/bash

#SBATCH --job-name=curl
#SBATCH --output=log/download_221014_A00643_0567_AHYMN3DSX3_P1713.log
#SBATCH --error=log/download_221014_A00643_0567_AHYMN3DSX3_P1713.log
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

# Global variables
link='https://filetransfer.mdc-berlin.de/?u=PVrmHjBp&p=dNNkQh7B'
files=$(cat docs/download/221014_A00643_0567_AHYMN3DSX3_P1713.txt)
out_dir="data/raw/tar/221014_A00643_0567_AHYMN3DSX3_P1713"

if [[ -d $out_dir ]]; then
  echo "$out_dir exists."
else
  mkdir $out_dir
fi

# Download files
for file in $files
do
  # Local variables
  url="$link&path=/$file"
  out_file="$out_dir/$file"

  echo "Download file from $url to $out_file"
  curl -C - -o $out_file $url
done

# Run md5sum on archives
archives=$(ls $out_dir/BCLS_LANE_*.tar)
for archive in $archives; do
  # Local variables
  out_file="$archive.md5loc"

  if [[ -e $out_file ]]; then
  echo "$out_file already exists and will be skipped."
  else
  echo "Running md5sum for $archive ..."
  md5sum $archive > $out_file
  fi

done
