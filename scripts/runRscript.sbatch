#!/bin/bash

#SBATCH --job-name=Rscript
#SBATCH -o SLURM_out/rscipt-%j.out
#SBATCH --time=4-00:00:00
#SBATCH -p normal,cbirdq
#SBATCH --nodes=1

to_run=$1
echo "${to_run}"

module load R/gcc/64/3.5.1


#https://stackoverflow.com/questions/3811345/how-to-pass-all-arguments-passed-to-my-bash-script-to-a-function-of-mine

Rscript "${to_run}" $2 #"${@:2}"
