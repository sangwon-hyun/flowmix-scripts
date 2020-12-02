#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --time=2-00:00:00
#SBATCH --export=none  # Ensures job gets a fresh login environment
#SBATCH --mail-user=robohyun66@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name=summarize
#SBATCH --output=./logs/summarize.%A_%a.log
#SBATCH --array=2-8
#SBATCH --exclusive

# not using this:   SBATCH --partition=scavenge

# # Set the session up to use R
# source /scratch/sangwonh/software/R/3.6.0/setup.sh
# # source /scratch/sangwonh/software/R/3.6.1/setup.sh

## New system:
module load gcc/8.3.0
module load openblas/0.3.7-oqk2
module load r/3.6.1-openblas-ydiq

ID="$SLURM_ARRAY_TASK_ID"
# eval "source /scratch/sangwonh/software/R/3.6.1/setup-${ID}.sh"


echo This is SLURM task $SLURM_ARRAY_TASK_ID

# Summarize the CV results by numclust (equal to SLURM array number)
Rscript summarize-numclust.R\
	nsim=30\
	numclust=$SLURM_ARRAY_TASK_ID\
