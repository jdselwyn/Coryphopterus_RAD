#!/bin/bash

#SBATCH --time=1-00:00:00
#SBATCH -p gpu,normal,cbirdq
#SBATCH --job-name=MitoID

module load samtools
module load bwa
module load mawk
module load ddocent
module load blast+/gcc7/2.9.0

#https://open.oregonstate.education/computationalbiology/chapter/command-line-blast/

REFERENCE="$(pwd)/$1" #Filepath to the reference mitochondrial genome fasta file ex) Reference_Sequence/bathygobius_cocosensis_complete_mitochondrion.fasta
BLASTDB="$(pwd)/$2" #Path to the BLAST DB ex) Reference_Sequence/CoryphopterusBlast
read_path="$(pwd)/$3" #Where the read files are stored (no "/" at the end) ex) mkREF_NovaSeq
output_path=$3 #Output Location ex) Mitochondrial_Mapping

readarray -t INDIVIDUALS < <(cut -d, -f2 ${output_path}/tmp_prefix_file) #Read in that individual indexing file
prefix=${INDIVIDUALS[${SLURM_ARRAY_TASK_ID}]}

echo Individual $prefix being mapped to $REFERENCE

#Need a way to do this if needed - or do this before the array
#samtools faidx ${REFERENCE}
#bwa index ${REFERENCE}

cd ${output_path}

mkdir -p bamfiles
mkdir -p fastqfiles
mkdir -p blast_out

bwa mem ${REFERENCE} ${read_path}/${prefix}*.r1.fq.gz ${read_path}/${prefix}*.r2.fq.gz -L 20,5 -t 20 -a -M -T 2 -A 1 -B 4 -O 6 -R "@RG\tID:${prefix}\tSM:${prefix}\tPL:Illumina" | mawk '!/\t[2-9].[SH].*/' | mawk '!/[2-9].[SH]\t/' | samtools view -@20 -q 1 -SbT ${REFERENCE} - > bamfiles/${prefix}.bam
samtools sort -@20 bamfiles/${prefix}.bam -o bamfiles/${prefix}.bam
samtools index bamfiles/${prefix}.bam

#Create fastq output
# samtools sort -@20 -n ${prefix}.bam -o ${prefix}_sort.bam
bedtools bamtofastq -i bamfiles/${prefix}.bam -fq fastqfiles/${prefix}.fastq
paste - - - - <fastqfiles/${prefix}.fastq | cut -f2 > fastqfiles/${prefix}.fasta

# Blast Results
touch blast_out/${prefix}.tsv

blastn \
  -query fastqfiles/${prefix}.fasta \
  -db ${BLASTDB} \
  -evalue 1e-6 \
  -outfmt 6 \
  -max_target_seqs 100 \
  -max_hsps 100 \
  -out blast_out/${prefix}.tsv \
  -num_threads ${SLURM_CPUS_ON_NODE}

rm tmp_prefix_file
