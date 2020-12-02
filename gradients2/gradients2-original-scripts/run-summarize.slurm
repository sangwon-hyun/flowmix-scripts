#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --time=7-00:00:00
#SBATCH --export=none  # Ensures job gets a fresh login environment
#SBATCH --mail-user=robohyun66@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name=summarize
#SBATCH --output=./logs/summarize.%A_%a.log
#SBATCH --array=80,81,82,83,84,85,86,89,87,88
#SBATCH --partition=scavenge
#SBATCH --exclusive

# 80-89

# not using this:   SBATCH --partition=scavenge

# Set the session up to use R
# source /staging/sh7/sangwonh/software/R/3.6.1/setup.sh

ID="$(($SLURM_ARRAY_TASK_ID - 79))"
eval "source /staging/sh7/sangwonh/software/R/3.6.1/setup-${ID}.sh"

echo This is SLURM task $SLURM_ARRAY_TASK_ID

# Main CV
Rscript summarize.R\
	datatype=$SLURM_ARRAY_TASK_ID\
	nsim=100\
