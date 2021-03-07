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

Rscript subsample-simulation.R\
	arraynum=$SLURM_ARRAY_TASK_ID\
	arraynum_max=$SLURM_ARRAY_TASK_MAX\
  	isim=$isim\
  	nsim=$nsim\
	subsample_size=$subsample_size\
	type=$type\
	summarize=$summarize\
	iboot=$iboot\

exit 0

#### Estimating convergence rate ####
nsim=50
for subsample_size in 20 40 80 100; do
    for isim in $(seq 1 1 $nsim); do
	jid1=$(sbatch --export=type=1,isim=$isim,iboot=$iboot,nsim=$nsim,subsample_size=$subsample_size,summarize=0 --array=1 run-subsample-simulation.slurm)
	jid2=$(sbatch --export=type=1,isim=$isim,iboot=$iboot,nsim=$nsim,subsample_size=$subsample_size,summarize=0 --array=1-20 --dependency=afterany:${jid1##* } run-subsample-simulation.slurm)
	jid3=$(sbatch --export=type=1,isim=$isim,iboot=$iboot,nsim=$nsim,subsample_size=$subsample_size,summarize=1 --dependency=afterany:${jid2##* } --array=1 run-subsample-simulation.slurm) ## experimental
    done
done


#### Run the actual simulations (iboot in 1:100)
nboot=50 ## NEED TO DO THIS IN CHUNKS.
nsim=100
subsample_size=50
for iboot in $(seq 1 1 $nboot); do
    for isim in $(seq 1 1 $nsim); do
	jid1=$(sbatch --export=type=2,isim=$isim,nsim=$nsim,subsample_size=$subsample_size,iboot=$iboot,summarize=0 --array=1 run-subsample-simulation.slurm)
	jid2=$(sbatch --export=type=2,isim=$isim,nsim=$nsim,subsample_size=$subsample_size,iboot=$iboot,summarize=0 --array=1-1 --dependency=afterany:${jid1##* } run-subsample-simulation.slurm)
	jid3=$(sbatch --export=type=2,isim=$isim,nsim=$nsim,subsample_size=$subsample_size,iboot=$iboot,summarize=1 --dependency=afterany:${jid2##* } --array=1 run-subsample-simulation.slurm)
    done
done

subsample_size=50
for iboot in 6 60; do
    for isim in $(seq 1 1 $nsim); do
	jid2=$(sbatch --export=type=2,isim=$isim,nsim=$nsim,subsample_size=$subsample_size,iboot=$iboot,summarize=0 --array=1-1 run-subsample-simulation.slurm)
	jid3=$(sbatch --export=type=2,isim=$isim,nsim=$nsim,subsample_size=$subsample_size,iboot=$iboot,summarize=1 --dependency=afterany:${jid2##* } --array=1 run-subsample-simulation.slurm)
    done
done

nboot=100
nsim=50
subsample_size=50
for iboot in $(seq 1 1 $nboot); do
    printf "bootstrap"
    printf $iboot
    printf "\n"
    printf "simulation "
    for isim in $(seq 1 1 $nsim); do
	printf $isim
	printf " "
	jid1=$(sbatch --export=type=3,isim=$isim,nsim=$nsim,subsample_size=$subsample_size,iboot=$iboot,summarize=0 --array=1 run-subsample-simulation.slurm)
	jid2=$(sbatch --export=type=3,isim=$isim,nsim=$nsim,subsample_size=$subsample_size,iboot=$iboot,summarize=0 --array=1-20 --dependency=afterany:${jid1##* } run-subsample-simulation.slurm)
	jid3=$(sbatch --export=type=3,isim=$isim,nsim=$nsim,subsample_size=$subsample_size,iboot=$iboot,summarize=1 --dependency=afterany:${jid2##* } --array=1 run-subsample-simulation.slurm)
    done
printf "\n"
done



# nsim=100
# subsample_size=50
# for iboot in $(seq 1 1 100); do
#     for isim in $(seq 1 1 $nsim); do
# 	sbatch --export=type=2,isim=$isim,nsim=$nsim,subsample_size=$subsample_size,iboot=$iboot,summarize=1 --array=1 run-subsample-simulation.slurm
#     done
# done



# exit
