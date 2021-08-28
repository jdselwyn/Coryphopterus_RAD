#!/bin/bash

# the hashtag and SBATCH are not commented out these are commands for the HPC
#SBATCH --time=1-00:00:00
#SBATCH -p gpu,normal,cbirdq
#SBATCH --nodes=1
#SBATCH --job-name=demultiArray

module load stacks

DIRECTORY=${1}

echo Hi
echo ${SLURM_ARRAY_TASK_ID}
echo ${DIRECTORY}


mapfile -t INDEXER < <(cat ${DIRECTORY}/fastq/index_array)
mapfile -t decode < <(cat ${DIRECTORY}/fastq/demulti_array)

echo ${INDEXER[${SLURM_ARRAY_TASK_ID}]}
echo ${decode[${SLURM_ARRAY_TASK_ID}]}

indexer=${INDEXER[${SLURM_ARRAY_TASK_ID}]}
DECODE=${decode[${SLURM_ARRAY_TASK_ID}]}
DECODE=$DIRECTORY/fastq/$DECODE

files=$(ls ${DIRECTORY}/fastq/*_[1-2].fq.gz | grep \_${indexer}\_)

SEQUENCES=${files[0]}
SEQ1=$(echo $SEQUENCES | grep -o '^[^ ]*')
SEQ2=$(echo $SEQUENCES | grep -o . | tac | tr -d "\n" | grep -o '^[^ ]*' | grep -o . | tac | tr -d "\n")


#create the directory where the demultiplexed sequences will go
mkdir -p ${DIRECTORY}/demultiplexed_seqs
mkdir -p ${DIRECTORY}/removed_seqs

mkdir -p ${DIRECTORY}/demultiplexed_seqs/${indexer}

process_radtags -1 ${SEQ1} -2 ${SEQ2} -i gzfastq -b ${DECODE} -r -e ecoRI -o ${DIRECTORY}/demultiplexed_seqs/${indexer} -D

# -1 designates the read 1 filel in a set of paired-end sequences. This needs to match the sequence file name.
# -2 designates the read 2 file in a set of paired-end sequences. This needs to match the sequence file name.
# -i sets the input file type. This should not be changed
# -b designates the barcode file that is used to separate the barcodes and name individual sequence files. This should be changed to match your barcode file name.
# -r tells the program to rescue barcodes that have no more than 2 mismatches. You do not need to adjust this parameter
# -e indicates the radtag that the program should search for. This is the enzyme associated with the read 1 sequences
# -o is the output path for the demultiplexed sequences. In this case it should be the demultiplexed seqs file
# -D tells the program to catch all of the discarded reads to a file

start="Barcode\tFilename"; end="^\s*$"; sed -n "/$start.*/,/$end/{/$start.*/b;/$end/b;p}" ${DIRECTORY}/demultiplexed_seqs/${indexer}/process_radtags.fastq.log | \
	sed -e 's/\t/,/g' > ${DIRECTORY}/demultiplexed_seqs/${indexer}_process_radtags.csv

# cleanup
rename 1.fq.gz F.fq.gz ${DIRECTORY}/demultiplexed_seqs/${indexer}/*
rename 2.fq.gz R.fq.gz ${DIRECTORY}/demultiplexed_seqs/${indexer}/*

rename radtags.fastq.log radtages_${indexer}.log ${DIRECTORY}/demultiplexed_seqs/${indexer}/process_radtags.fastq.log
mv ${DIRECTORY}/demultiplexed_seqs/${indexer}/* ${DIRECTORY}/demultiplexed_seqs/
rmdir ${DIRECTORY}/demultiplexed_seqs/${indexer}

#mkdir -p $DIRECTORY/demultiplexed_seqs/removed_seqs
mv ${DIRECTORY}/demultiplexed_seqs/*.rem.* ${DIRECTORY}/removed_seqs
