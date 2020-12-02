#!/bin/bash

nsim=100
# 80 81 82 83 84 85 86 87 88
## run the following:
##81 82 84 85 87 88
for datatype in 80 81 82 83 84 85 86 87 88 89
do  echo "running data type:"
    echo $datatype
    # sbatch -J "1-${datatype}-2" run-blockcv-1-8x-2.slurm $datatype
    sbatch -J "2-${datatype}-2" run-sim.slurm $datatype $nsim
done
