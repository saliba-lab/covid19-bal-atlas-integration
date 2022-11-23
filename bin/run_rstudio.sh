#!/bin/bash

# Define variables
base_dir=$(pwd)
env=covid19-bal-atlas
conda_path=/home/$USER/miniconda3

# Activate local environment
export PATH=$conda_path/envs/$env/bin:$PATH
export LD_LIBRARY_PATH=$conda_path/lib

# Run rstudio
rstudio &

