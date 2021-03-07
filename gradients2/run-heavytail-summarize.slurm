#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --export=none  # Ensures job gets a fresh login environment
#SBATCH --mail-user=robohyun66@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name=heavytail-summarize
#SBATCH --output=./logs/heavytail-summarize.%A_%a.log
#SBATCH --error=./logs/heavytail-summarize.%A_%a.err

# Set the session up to use R
module load gcc/8.3.0
module load openblas/0.3.8
module load r/3.6.3

echo This is SLURM task $SLURM_ARRAY_TASK_ID

Rscript heavytail-summarize.R\
	arraynum=$SLURM_ARRAY_TASK_ID\
	arraynum_max=$SLURM_ARRAY_TASK_MAX\
  	isim=$isim\
  	nsim=$nsim\
	ialpha=$ialpha\
	DF=$DF\
	noisetype_ii=$noisetype_ii\

exit

## Gaussian
nsim=200
ialpha=10000
noisetype_ii=1
DF=-10000
for isim in $(seq 1 1 $nsim); do
    sbatch --export=isim=$isim,nsim=$nsim,ialpha=$ialpha,DF=$DF,noisetype_ii=$noisetype_ii --array=1 run-heavytail-summarize.slurm
done

## Use this to run the heavy tailed experiments
nsim=200
ialpha=10000
noisetype_ii=2
for DF in 2 5 10 20 40 100 ; do
    for isim in $(seq 1 1 $nsim); do
	sbatch --export=isim=$isim,nsim=$nsim,ialpha=$ialpha,DF=$DF,noisetype_ii=$noisetype_ii --array=1 run-heavytail-summarize.slurm
    done
done


## Use this to run the skewed data experiments
nsim=200
noisetype_ii=1
DF=-10000
for ialpha in 1 2 3 4 5 ; do
    for isim in $(seq 1 1 $nsim); do
	sbatch --export=isim=$isim,nsim=$nsim,ialpha=$ialpha,DF=$DF,noisetype_ii=$noisetype_ii --array=1 run-heavytail-summarize.slurm
    done
done
