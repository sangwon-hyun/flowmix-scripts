#!/bin/bash

## 3 5 7
##2 4 6 8 9 10
# 5 7 4 6 9 2 8 10
nsim=30
for numclust in 8 7
do  echo "running data type:"
    echo $numclust
    # sbatch  -J "2-9-${numclust}" run-blockcv-2-9-x-discovery.slurm $numclust $nsim
    sbatch  -J "2-9-${numclust}" run-blockcv-2-9-x-discovery-mix.slurm $numclust $nsim
    # sbatch  -J "2-9-${numclust}" run-blockcv-2-9-x.slurm $numclust $nsim
done
