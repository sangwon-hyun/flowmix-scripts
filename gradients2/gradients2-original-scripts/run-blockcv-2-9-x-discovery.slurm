#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --time=2-00:00:00
#SBATCH --export=none  # Ensures job gets a fresh login environment
#SBATCH --mail-user=robohyun66@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=./logs/blockcv-2-9-x.%A_%a.log
#SBATCH --array=1-20
#SBATCH --exclusive

# Set the session up to use R
# source /staging/sh7/sangwonh/software/R/3.6.1/setup.sh
# eval "source /staging/sh7/sangwonh/software/R/3.6.1/setup-${SLURM_ARRAY_TASK_ID}.sh"

echo Set up starting

# Original system:
# eval "source /scratch/sangwonh/software/R/3.6.0/setup.sh"

## New system:
module load gcc/8.3.0
module load openblas/0.3.7-oqk2
module load r/3.6.1-openblas-ydiq

# module load r/3.6.3

echo Set up ended

echo This is SLURM task $SLURM_ARRAY_TASK_ID

# Main CV
Rscript blockcv-numclust.R\
	arraynum=$SLURM_ARRAY_TASK_ID\
	arraynum_max=20\
	blocktype=2\
        datatype=9\
	numclust=$1\
	cv_gridsize=7\
	nfold=5\
	nrep=10\
	maxres_once=0\
	sim=1\
	nsim=$2\
