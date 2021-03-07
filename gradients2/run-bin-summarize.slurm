#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --export=none  # Ensures job gets a fresh login environment
#SBATCH --mail-user=robohyun66@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name=bin-summarize
#SBATCH --output=./logs/bin-summarize.%A_%a.log
#SBATCH --error=./logs/bin-summarize.%A_%a.err

# Set the session up to use R
module load gcc/8.3.0
module load openblas/0.3.8
module load r/3.6.3

echo This is SLURM task $SLURM_ARRAY_TASK_ID

Rscript bin-summarize.R\
	arraynum=$SLURM_ARRAY_TASK_ID\
	arraynum_max=$SLURM_ARRAY_TASK_MAX\
  	isim=$isim\
  	nsim=$nsim\
	B=$bin_size\

exit

## Use this to run the main experiments
nsim=100
for bin_size in 5 10 20 30 40 50 100 10000; do
#     # Estimate the maximum lambda value once.
    for isim in $(seq 1 1 $nsim); do
	sbatch --export=isim=$isim,nsim=$nsim,bin_size=$bin_size --array=1-200 run-bin-summarize.slurm
    done
done
