#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --time=7-00:00:00
#SBATCH --export=none  # Ensures job gets a fresh login environment
#SBATCH --mail-user=robohyun66@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=./logs/blockcv-2-9-x-fi.%A_%a.log
#SBATCH --array=1
#SBATCH --partition=scavenge
#SBATCH --exclusive
#SBATCH --constraint=IB

# Set the session up to use R
# source /staging/sh7/sangwonh/software/R/3.6.1/setup.sh
# eval "source /staging/sh7/sangwonh/software/R/3.6.1/setup-${SLURM_ARRAY_TASK_ID}.sh"

echo Set up starting

# eval "source /scratch/sangwonh/software/R/3.6.1/setup-${SLURM_ARRAY_TASK_ID}.sh"
# eval "source /scratch/sangwonh/software/R/3.6.1/setup.sh"
eval "source /scratch/sangwonh/software/R/3.6.0/setup.sh"

echo Set up ended

echo This is SLURM task $SLURM_ARRAY_TASK_ID

# Main CV
Rscript blockcv.R\
	arraynum=$SLURM_ARRAY_TASK_ID\
	arraynum_max=50\
	blocktype=2\
        datatype=9\
	numclust=$1\
	cv_gridsize=7\
	nfold=5\
	nrep=5\
	random_order=0\
	maxres_once=0\
	fitonly=0\
	high_range_beta=0\
	high_range=0\
	sim=1\
	nsim=$2\

# Fit only
Rscript blockcv.R\
	arraynum=$SLURM_ARRAY_TASK_ID\
	arraynum_max=50\
	blocktype=2\
	datatype=9\
	numclust=$1\
	cv_gridsize=7\
	nfold=5\
	nrep=5\
	random_order=0\
	maxres_once=0\
	fitonly=1\
	high_range_beta=0\
	high_range=0\
	sim=1\
	nsim=$2\
