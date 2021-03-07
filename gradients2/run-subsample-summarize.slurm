#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --export=none  # Ensures job gets a fresh login environment
#SBATCH --mail-user=robohyun66@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name=subsample-summarize
#SBATCH --output=./logs/subsample-summarize.%A_%a.log
#SBATCH --error=./logs/subsample-summarize.%A_%a.err

# Set the session up to use R
module load gcc/8.3.0
module load openblas/0.3.8
module load r/3.6.3

echo This is SLURM task $SLURM_ARRAY_TASK_ID

Rscript subsample-summarize.R\
	arraynum=$SLURM_ARRAY_TASK_ID\
	arraynum_max=$SLURM_ARRAY_TASK_MAX\
  	isim=$isim\
  	nsim=$nsim\
	subsample_size=$subsample_size\

exit

## Use this to run the main experiments
nsim=100
for subsample_size in 30 60 100 148 200 250; do
#     # Estimate the maximum lambda value once.
    for isim in $(seq 1 1 $nsim); do
	jid1=$(sbatch --export=isim=$isim,nsim=$nsim,subsample_size=$subsample_size --array=1 run-subsample.slurm)
	sbatch --export=isim=$isim,nsim=$nsim,subsample_size=$subsample_size --array=1-20 --dependency=afterany:${jid1##* } run-subsample.slurm
    done
done
