#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --export=none  # Ensures job gets a fresh login environment
#SBATCH --mail-user=robohyun66@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name=heavytail
#SBATCH --output=./logs/heavytail.%A_%a.log
#SBATCH --error=./logs/heavytail.%A_%a.err

# Set the session up to use R
module load gcc/8.3.0
module load openblas/0.3.8
module load r/3.6.3

echo This is SLURM task $SLURM_ARRAY_TASK_ID

Rscript heavytail.R\
	arraynum=$SLURM_ARRAY_TASK_ID\
	arraynum_max=$SLURM_ARRAY_TASK_MAX\
  	isim=$isim\
  	nsim=$nsim\
	ialpha=$ialpha\
	DF=$DF\
	noisetype_ii=$noisetype_ii\
	run_summarize=$run_summarize\

exit

## Gaussian
nsim=200
ialpha=10000
noisetype_ii=1
DF=-10000
for isim in $(seq 1 1 $nsim); do
    # Estimate the maximum lambda value once.
    jid1=$(sbatch --export=isim=$isim,nsim=$nsim,ialpha=$ialpha,DF=$DF,noisetype_ii=$noisetype_ii,run_summarize=0 --array=1 run-heavytail.slurm)
    # Run the actual experiments
    jid2=$(sbatch --export=isim=$isim,nsim=$nsim,ialpha=$ialpha,DF=$DF,noisetype_ii=$noisetype_ii,run_summarize=0 --array=1-10 --dependency=afterany:${jid1##* } run-heavytail.slurm)
    # Summarize
    jid3=$(sbatch --export=isim=$isim,nsim=$nsim,ialpha=$ialpha,DF=$DF,noisetype_ii=$noisetype_ii,run_summarize=1 --array=1 --dependency=afterany:${jid2##* } run-heavytail.slurm)
done

## Use this to run the heavy tailed experiments
nsim=100
ialpha=10000
noisetype_ii=2
for DF in 3 5 10 20 40 100 ; do
# for DF in 3 5 10 20 40 100 ; do
    for isim in $(seq 1 1 $nsim); do
	# Estimate the maximum lambda value once.
	jid1=$(sbatch --export=isim=$isim,nsim=$nsim,ialpha=$ialpha,DF=$DF,noisetype_ii=$noisetype_ii,run_summarize=0 --array=1 run-heavytail.slurm)
	# Run the actual experiments
	jid2=$(sbatch --export=isim=$isim,nsim=$nsim,ialpha=$ialpha,DF=$DF,noisetype_ii=$noisetype_ii,run_summarize=0 --array=1-2 --dependency=afterany:${jid1##* } run-heavytail.slurm)
	# Summarize
	jid3=$(sbatch --export=isim=$isim,nsim=$nsim,ialpha=$ialpha,DF=$DF,noisetype_ii=$noisetype_ii,run_summarize=1 --array=1 --dependency=afterany:${jid2##* } run-heavytail.slurm)
    done
done


## Use this to run the skewed data experiments
nsim=100
noisetype_ii=3
DF=10000
for ialpha in 1 2 3 4 5 ; do
    for isim in $(seq 1 1 $nsim); do
	# Estimate the maximum lambda value once.
	jid1=$(sbatch --export=isim=$isim,nsim=$nsim,ialpha=$ialpha,DF=$DF,noisetype_ii=$noisetype_ii,run_summarize=0 --array=1 run-heavytail.slurm)
	# Run the actual experiments
	jid2=$(sbatch --export=isim=$isim,nsim=$nsim,ialpha=$ialpha,DF=$DF,noisetype_ii=$noisetype_ii,run_summarize=0 --array=1 --dependency=afterany:${jid1##* } run-heavytail.slurm)
	# Summarize
	jid3=$(sbatch --export=isim=$isim,nsim=$nsim,ialpha=$ialpha,DF=$DF,noisetype_ii=$noisetype_ii,run_summarize=1 --array=1 --dependency=afterany:${jid2##* } run-heavytail.slurm)
    done
done
