#!/bin/bash

#SBATCH --time=1-00:00:00
#SBATCH -p gpu,normal,cbirdq
#SBATCH --job-name=mito_map

module load samtools
module load bwa
module load mawk
module load ddocent
module load qualimap

REFERENCE=$1 #First argument should be the filepath to the reference mitochondrial genome fasta file ex) /work/hobi/jselwyn/MitoMap/lottia_digitalis_mito_sequence.fasta
read_path=$2 #Second argument should be the path to where the read files are stored (no "/" at the end)
output_path=$3 #Third argument is the final location to put the output files

readarray -t INDIVIDUALS < <(cut -d, -f2 ${output_path}/tmp_prefix_file) #Read in that individual indexing file
prefix=${INDIVIDUALS[${SLURM_ARRAY_TASK_ID}]}

echo Individual $prefix being mapped to $REFERENCE

#Need a way to do this if needed - or do this before the array
#samtools faidx ${REFERENCE}
#bwa index ${REFERENCE}

cd ${output_path}

mkdir -p bamfiles
mkdir -p fastqfiles

bwa mem ${REFERENCE} ${read_path}/${prefix}.r1.fq.gz ${read_path}/${prefix}.r2.fq.gz -L 20,5 -t 20 -a -M -T 2 -A 1 -B 4 -O 6 -R "@RG\tID:${prefix}\tSM:${prefix}\tPL:Illumina" | mawk '!/\t[2-9].[SH].*/' | mawk '!/[2-9].[SH]\t/' | samtools view -@20 -q 1 -SbT ${REFERENCE} - > bamfiles/${prefix}.bam
samtools sort -@20 bamfiles/${prefix}.bam -o bamfiles/${prefix}.bam
samtools index bamfiles/${prefix}.bam

#Create fastq output
# samtools sort -@20 -n ${prefix}.bam -o ${prefix}_sort.bam
bedtools bamtofastq -i bamfiles/${prefix}.bam -fq fastqfiles/${prefix}.fastq
