#!/bin/bash

# Set up variables
base_dir=$(pwd)
env=covid19-bal-atlas
conda_path=/home/$USER/miniconda3

# Activate local environment
export PATH=$conda_path/envs/$env/bin:$PATH

# Run Spyder
spyder &
