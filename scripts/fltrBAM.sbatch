#!/bin/bash

#SBATCH --job-name=fltrBAM
#SBATCH -o SLURM_out/fltrBAM-%j.out
#SBATCH --time=48:00:00
#SBATCH -p normal
#SBATCH --nodes=1

module load ddocent

directory=${1}
configFile=${2}

dDocent="$(pwd)/scripts/dDocentHPC/dDocentHPC.bash"
cp ${configFile} ${directory}/fltrBAM.config

#this will use dDocent fltrBAM to filter the BAM files
cd ${directory}
bash ${dDocent} fltrBAM fltrBAM.config
