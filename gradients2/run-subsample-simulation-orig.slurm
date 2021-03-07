#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --export=none  # Ensures job gets a fresh login environment
#SBATCH --mail-user=robohyun66@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name=subsample-sim
#SBATCH --output=./logs/subsample-sim.%A_%a.log
#SBATCH --error=./logs/subsample-sim.%A_%a.err

# Set the session up to use R
module load gcc/8.3.0
module load openblas/0.3.8

module load r/3.6.3

echo This is SLURM task $SLURM_ARRAY_TASK_ID

Rscript subsample-simulation-orig.R\
	arraynum=$SLURM_ARRAY_TASK_ID\
	arraynum_max=$SLURM_ARRAY_TASK_MAX\
	summarize=$summarize\
	iboot=$iboot\

exit 0


### Original (no subsampling, full-data) estimation
nboot=100
for iboot in $(seq 2 1 $nboot); do
	jid1=$(sbatch --export=iboot=$iboot,summarize=0 --array=1 run-subsample-simulation-orig.slurm)
	jid2=$(sbatch --export=iboot=$iboot,summarize=0 --array=1-20 --dependency=afterany:${jid1##* } run-subsample-simulation-orig.slurm)
	jid3=$(sbatch --export=iboot=$iboot,summarize=1 --dependency=afterany:${jid2##* } --array=1 run-subsample-simulation-orig.slurm) ## experimental
    done
done

