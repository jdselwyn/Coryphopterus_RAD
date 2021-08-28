#!/bin/bash

#SBATCH --job-name=mkREF
#SBATCH -o SLURM_out/mkREF-%j.out
#SBATCH --time=4-00:00:00
#SBATCH -p normal,cbirdq
#SBATCH --nodes=1

module load ddocent/2.7.8

directory=${1}
configFile=${2}

dDocent="$(pwd)/scripts/dDocentHPC/dDocentHPC.bash"

cp ${configFile} ${directory}/mkREF.config

cd ${directory}

bash ${dDocent} mkREF mkREF.config
