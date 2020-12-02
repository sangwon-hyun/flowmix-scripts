#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --export=none  # Ensures job gets a fresh login environment
#SBATCH --mail-user=robohyun66@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name=gradients1
#SBATCH --output=./logs/gradients1.%A_%a.log
#SBATCH --array=1-1200

# not using this: SBATCH --exclusive
# not using this:   SBATCH --partition=scavenge

# Set the session up to use R
module load gcc/8.3.0
module load openblas/0.3.8
module load r/3.6.3

echo This is SLURM task $SLURM_ARRAY_TASK_ID

# Main CV
Rscript gradients1.R\
	arraynum=$SLURM_ARRAY_TASK_ID\
	arraynum_max=1200\
	maxres=0\
	fitonly=0\
	numclust=7\
	blocksize=5\
	nrep=20\

Rscript gradients1.R\
	arraynum=$SLURM_ARRAY_TASK_ID\
	arraynum_max=1200\
	maxres=0\
	fitonly=1\
	numclust=7\
	blocksize=5\
	nrep=20\


# # Main CV
# Rscript gradients1.R\
# 	arraynum=$SLURM_ARRAY_TASK_ID\
# 	arraynum_max=1200\
# 	maxres=0\
# 	fitonly=1\
