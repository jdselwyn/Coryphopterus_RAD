#!/bin/bash

#SBATCH --time=4-00:00:00
#SBATCH -p normal
#SBATCH --job-name=refIndex
#SBATHC --output=SLURM_Out/index_ref_%j.out

#Ex usage
#sbatch index_ref.slurm $WORK/Genetic_Data/NovaSeq/Reference_CO1/bathygobius_cocosensis_complete_mitochondrion.fasta

module load samtools
module load bwa

REFERENCE=$1

samtools faidx ${REFERENCE}
bwa index ${REFERENCE}
