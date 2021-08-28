#!/bin/bash

#SBATCH --job-name=repr_fqscrn
#SBATCH -o SLURM_out/repair_fastqscreen_%A.out
#SBATCH -p normal,cbirdq

module load java/1.8
module load bbmap/38.90
module load parallel

INDIR=$1
OUTDIR=$2
THREADS=${SLURM_CPUS_ON_NODE}
TEMPDIR=/tmp
EXTENSION=.tagged_filter.fastq.gz
#EXTENSION2=.tagged_filter.fastq

mkdir $OUTDIR $OUTDIR/orphans

# unzip the fq.gz files
#ls $INDIR/*$EXTENSION | \
#	sed -e 's/.*\///g' | \
#	parallel --no-notice -j $THREADS \
#	gunzip $INDIR/{}

# repair the paired end files
ls $INDIR/*$EXTENSION | \
	sed -e "s/r[12]$EXTENSION//" -e "s/orph$EXTENSION//" -e 's/.*\///g' | \
	uniq | \
	parallel --no-notice -j $THREADS \
		repair.sh \
		in1=$INDIR/{}r1$EXTENSION \
		in2=$INDIR/{}r2$EXTENSION \
		out1=$OUTDIR/{}repr.r1.TMP.fq.gz \
		out2=$OUTDIR/{}repr.r2.TMP.fq.gz \
		outs=$OUTDIR/orphans/{}repr.orph.fq.gz \
		repair

#Fix issue if there is no space between "/[12]" and "#" seperating which read it is and fltrscrn
ls $OUTDIR/*TMP.fq.gz | \
  sed -e "s/TMP.fq.gz//" -e 's/.*\///g' | \
  parallel --gnu --no-notice -j $THREADS \
    "gunzip -c $OUTDIR/{}TMP.fq.gz | \
    sed -e '/^@.*\/[12]#/{s!#! #!g}' | \
    gzip > $OUTDIR/{}fq.gz"

rm $OUTDIR/*TMP.fq.gz


#ls $OUTDIR/*fq | \
#	parallel --no-notice -j $THREADS \
#		gzip {}

# concat orphans created by repair with the orphans from previous operations on pe files and save to outdir
# ls $INDIR/*$EXTENSION2 | \
# sed -e "s/r[12]$EXTENSION2//" -e "s/orph$EXTENSION2//" -e 's/.*\///g' | \
# uniq | \
# parallel --no-notice -j $THREADS \
        # "cat $INDIR/{}orph$EXTENSION2 $OUTDIR/{}repr_orph.fq > $OUTDIR/{}orph.fq"

# # remove orphans that are not concatted in outdir
# ls $OUTDIR/*repr_orph.fq | \
# crun parallel --no-notice -j40 \
        # rm {}
