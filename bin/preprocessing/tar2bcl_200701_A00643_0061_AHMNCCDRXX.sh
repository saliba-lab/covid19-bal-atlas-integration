#!/bin/bash

# set parameters for SGE submission
#$ -N count
#$ -l arch=linux-x64
#$ -pe multislot 2
#$ -b n
#$ -q all.q
#$ -o log/tar2bcl.log
#$ -e log/tar2bcl.log
#$ -cwd

# Print hostname & date
hostname -f
date

# Flowcell AHMNCCDRXX
from="data/raw/tar/200701_A00643_0061_AHMNCCDRXX.tar"
to="data/raw/bcl"
tar -xvf $from --directory $to
