#!/bin/bash

#SBATCH --job-name=mkVCF
#SBATCH -o SLURM_out/mkVCF-%j.out
#SBATCH --time=4-00:00:00
#SBATCH -p normal,cbirdq
#SBATCH --nodes=1

module load ddocent/2.7.8

directory=${1}
configFile=${2}
refGenome=${3}

dDocent="$(pwd)/scripts/dDocentHPC/dDocentHPC.bash"

#Config cutoffs in mkREF section must match the reference genome being used
cp ${configFile} ${directory}/mkVCF.config
cp ${refGenome} ${directory}

#this will use dDocent mkBAM to map the reads in the mkBAM directory to the reference*fasta created by
#dDocent mkREF.  We have to change directories, copy the reference*fasta to mkBAM, then run dDocent mkBAM
#
cd ${directory}
bash ${dDocent} mkVCF mkVCF.config
