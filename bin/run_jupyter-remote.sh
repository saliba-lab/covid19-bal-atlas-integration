#!/bin/bash
#SBATCH --job-name=jupyter_remote
#SBATCH --output=log/jupyter_remote.log
#SBATCH --error=log/jupyter_remote.log
#SBATCH --ntasks=2
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
# if GPUs needed: #SBATCH --partition=gpu
# if GPUs needed: #SBATCH --gres=gpu:t4:2
#SBATCH --time=48:00:00
#SBATCH --mem=256G

# Conda initialize #########################################

__conda_setup="$('/home/$USER/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/$USER/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/$USER/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/$USER/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup

# Variables ################################################
env="covid19-bal-atlas"
port="8080"

# Main program #############################################
unset PYTHONPATH
conda activate $env
jupyter notebook --no-browser --port=$port
