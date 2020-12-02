#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --export=none  # Ensures job gets a fresh login environment
#SBATCH --mail-user=robohyun66@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name=gradients2-sim-noisy-X
#SBATCH --output=./logs/gradients2-sim-noisy-X.%A_%a.log
#SBATCH --array=1-1200

# Set the session up to use R
module load gcc/8.3.0
module load openblas/0.3.8
module load r/3.6.3

echo This is SLURM task $SLURM_ARRAY_TASK_ID

# Main CV
Rscript gradients2.R\
	arraynum=$SLURM_ARRAY_TASK_ID\
	arraynum_max=1200\
	experiment_type=4\
	numclust=5\
	cv_gridsize=10\
	nfold=5\
	nrep=10\
	random_order=0\
	maxres_once=0\
	fitonly=0\
	sim=1\
	nsim=100\


# Fit only
Rscript gradients2.R\
	arraynum=$SLURM_ARRAY_TASK_ID\
	arraynum_max=1200\
	experiment_type=4\
	numclust=5\
	cv_gridsize=10\
	nfold=5\
	nrep=10\
	random_order=0\
	maxres_once=0\
	fitonly=1\
	sim=1\
	nsim=100\
