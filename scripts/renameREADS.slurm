#!/bin/bash

#SBATCH --job-name=renameReads
#SBATCH -o /work/hobi/jselwyn/Coryphopterus_RAD/testNovaSeq_Assembly/renameReads_%j.out
#SBATCH -p normal,cbirdq

module load parallel

cd /work/hobi/jselwyn/Coryphopterus_RAD/testNovaSeq_Assembly

THREADS=${SLURM_CPUS_ON_NODE}

ls *clmp.fp2.repr*fq.gz | \
  sed -e "s/fq.gz//" -e 's/.*\///g' | \
  parallel --gnu --no-notice -j $THREADS \
    "gunzip -c {}fq.gz | \
    sed -e '/^@.*\/[12]/{s! .*!!g}' | \
    gzip > {}renamed_fq.gz"
