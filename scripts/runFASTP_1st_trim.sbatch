#!/bin/bash

# this script will do all trimming, except 5'
# no merging of overlapping reads
# this is first step in prepping reads for de novo assembly

#SBATCH --job-name=fp1nss
#SBATCH -o SLURM_out/fastp_1st_noSizeSelect-%j.out
#SBATCH --time=4-00:00:00
#SBATCH -p normal,cbirdq

module load fastp/gcc7/0.20.1
module load parallel
module load python3/gcc7/3.7.4 #multiqc

INDIR=${1}
OUTDIR=${2}
MINLEN=${3}  # delete reads with fewer nt after trimming
FQPATTERN=*.fq.gz  #determines files to be trimmed
EXTPATTERN=[RF]\.fq\.gz  # pattern match to fq extensions
FWDEXT=F.fq.gz
REVEXT=R.fq.gz
THREADS=${SLURM_CPUS_ON_NODE}

mkdir $OUTDIR $OUTDIR/failed

#Trim
ls $INDIR/$FQPATTERN | \
  sed -e "s/$EXTPATTERN//" -e 's/.*\///g' | \
  uniq | \
  parallel --no-notice -j $THREADS \
  fastp \
    --in1 $INDIR/{}$FWDEXT \
    --in2 $INDIR/{}$REVEXT \
    --out1 $OUTDIR/{}r1.fq.gz \
    --out2 $OUTDIR/{}r2.fq.gz \
    --unpaired1 $OUTDIR/{}unprd.fq.gz \
    --unpaired2 $OUTDIR/{}unprd.fq.gz \
    --failed_out $OUTDIR/failed/{}fail.fq.gz \
    -h $OUTDIR/{}fastp.html \
    -j $OUTDIR/{}fastp.json \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 40 \
    --length_required $MINLEN \
    --low_complexity_filter \
    --complexity_threshold 30 \
    --detect_adapter_for_pe \
    --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --cut_tail \
    --cut_tail_window_size 1 \
    --cut_tail_mean_quality 20 \
    --trim_poly_g \
    --poly_g_min_len 10 \
    --trim_poly_x \
    --report_title "fastp_1"

mv $OUTDIR/*unprd* $OUTDIR/failed

#add multiqc line when it's updated
multiqc $OUTDIR -n $OUTDIR/1st_fastp_report
