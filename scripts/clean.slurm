#!/bin/bash

#SBATCH --job-name=compress
#SBATCH -p normal,cbirdq
#SBATCH --nodes=1
#SBATCH -t 4-00:00:00
#SBATCH -o SLURM_out/compress_%j.out

#### compress files ####
inFolder=$1
outFolder=$2

tar -cjvf ${outFolder}.tar.bz ${inFolder}
#sbatch scripts/clean.slurm ../prj_garcia_nudibranchs /work/hobi/webshare/prj_garcia_nudibranchs/prj_garcia_nudibranchs_V1
