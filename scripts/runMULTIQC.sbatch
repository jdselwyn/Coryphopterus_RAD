#!/bin/bash

#SBATCH --job-name=mltqc
#SBATCH -o SLURM_out/multiqc-%j.out
#SBATCH --time=24:00:00
#SBATCH -p normal,gpu,cbirdq

module load python3/gcc7/3.7.4 #multiqc

OUTDIR=$1

#add multiqc line when it's updated
multiqc $OUTDIR -n $OUTDIR/multiqc_report
