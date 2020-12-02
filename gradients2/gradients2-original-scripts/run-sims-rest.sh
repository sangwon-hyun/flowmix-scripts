#!/bin/bash
for datatype in 80 81 82 83 84 85 86 87 88 89
do  echo "running data type:"
    echo $datatype
    sbatch -J "2-${datatype}-2" run-blockcv-2-8x-2.slurm $datatype
done


# 83 84 85 86 87 88 89
