#!/bin/bash

# sbatch -J "sim-maxres" run-sim-maxres.slurm

# # This script assumes that the script has maxres=TRUE


# for datatype in 80 81 82 83 84 85 86 87 88 89
# # for datatype in 80 83 86 89
# do  echo "running data type:"
#     echo $datatype
#     # sbatch -J "1-${datatype}-2" run-blockcv-1-8x-2-maxres.slurm $datatype
#     sbatch -J "2-${datatype}-2" run-blockcv-2-8x-2-maxres.slurm $datatype
# done

# for numclust in 8
for numclust in 2 3 4 5 6 7 8
do  echo "running numclust:"
    echo $numclust
    sbatch -J "2-9-${numclust}" run-blockcv-2-9-x-maxres.slurm $numclust
done
