#!/bin/bash

# this script runs fastq_screen and outputs a report AND filtered data
# `mapdamage_log.xlsx` will create the arguments which should be copied into a config file

#SBATCH --job-name=fqscr_array

module load parallel/gcc7/20190222
module load fastq_screen/0.14.1
module load bowtie2/gcc7/2.3.5.1

FQPATTERN=$1       					#PATTERN="*.fq.gz"
INDIR=$2       					#INDIR=/work/hobi/cbird/pire_shotgunWGAvNoWGA/4sumstats_fq_fstp
TMPDIR=/tmp/fastqscrn_${SLURM_JOB_USER}_${SLURM_JOBID}_${SLURM_ARRAY_TASK_ID}      							#TMPDIR=/tmp/fastqscrn
OUTDIR=$3        	#OUTDIR=/work/hobi/cbird/pire_shotgunWGAvNoWGA/4sumstats_fq_fstp_fqscrn
CONFFILE=/work/hobi/egarcia111/reference_databases/fastq_screen_nofish_20210424.conf      					#CONFFILE=/work/hobi/egarcia111/reference_databases/fastq_screen_nofish_20210424.conf
SUBSET=0
FILTER=000000000000
ALIGNER=bowtie2  #bwa or bowtie2
THREADS=${SLURM_CPUS_ON_NODE}     # ram usage can be high, >100gb, so don't increase this beyond PARALLELISM * 120gb < ram avail

echo INDIR=$INDIR
echo TMPDIR=$TMPDIR
echo OUTDIR=$OUTDIR
echo CONFFILE=$CONFFILE
echo PATTERN=$FQPATTERN
echo PARALLELISM=$PARALLELISM
echo THREADS=$THREADS
echo SUBSET=$SUBSET
echo FILTER=$FILTER


all_samples=( $(ls $INDIR/$FQPATTERN) )
INFILES=${all_samples[${SLURM_ARRAY_TASK_ID}]}
echo INFILES:
echo ${INFILES}

rm -rf $TMPDIR
mkdir $TMPDIR

# if data is being subset, then don't output fqgz files
if [ $SUBSET -gt 0 ]; then
	echo ""; echo `date` RUNNING FASTQ_SCREEN TO GENERATE REPORTS, SUBSET=$SUBSET
# else output fqgz files
else
	echo ""; echo `date` RUNNING FASTQ_SCREEN TO FILTER FQ FILES AND GENERATE REPORTS, NO SUBSETTING
fi

fastq_screen \
	--aligner $ALIGNER \
	--conf $CONFFILE \
	--threads $THREADS \
	--tag \
	--force \
	--filter $FILTER \
	--subset $SUBSET \
	--outdir $TMPDIR \
	${INFILES}

ls $TMPDIR/* | \
	parallel --no-notice -j ${THREADS} \
		cp {} $OUTDIR

rm -rf $TMPDIR
