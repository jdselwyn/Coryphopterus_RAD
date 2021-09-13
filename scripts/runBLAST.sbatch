#!/bin/bash

#SBATCH --time=1-00:00:00
#SBATCH -p gpu,normal,cbirdq
#SBATCH --job-name=BLAST

#https://open.oregonstate.education/computationalbiology/chapter/command-line-blast/

module load blast+/gcc7/2.9.0

BLASTDB="$(pwd)/$1" #Reference_Sequence/CoryphopterusBlast
read_path="$(pwd)/${2}" #Mitochondrial_Mapping

readarray -t INDIVIDUALS < <(cut -d, -f2 ${read_path}/tmp_prefix_file) #Read in that individual indexing file
prefix=${INDIVIDUALS[${SLURM_ARRAY_TASK_ID}]}

echo Individual $prefix being BLASTED against $BLASTDB

mkdir -p ${read_path}/blast_out

paste - - - - <${read_path}/fastqfiles/${prefix}.fastq | cut -f2 > ${read_path}/fastqfiles/${prefix}.fasta

touch ${read_path}/blast_out/${prefix}.tsv

blastn \
  -query ${read_path}/fastqfiles/${prefix}.fasta \
  -db ${BLASTDB} \
  -evalue 1e-6 \
  -outfmt 6 \
  -max_target_seqs 100 \
  -max_hsps 100 \
  -out ${read_path}/blast_out/${prefix}.tsv \
  -num_threads ${SLURM_CPUS_ON_NODE}
