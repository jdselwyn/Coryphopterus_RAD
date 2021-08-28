#!/bin/bash

# this script will do all trimming, except 5'
# no merging of overlapping reads
# this is first step in prepping reads for de novo assembly

#SBATCH --job-name=fp2
#SBATCH -o SLURM_out/fastp_2nd-%j.out
#SBATCH --time=96:00:00
#SBATCH -p normal,cbirdq

module load fastp/gcc7/0.20.1
module load parallel
module load python3/gcc7/3.7.4 #multiqc

INDIR=${1}
OUTDIR=${2}
MINLEN=${3}  # delete reads with fewer nt after trimming
FQPATTERN=*.fq.gz  #determines files to be trimmed
EXTPATTERN=r[12]\.fq\.gz  # pattern match to fq extensions
FWDEXT=r1.fq.gz
REVEXT=r2.fq.gz
THREADS=${SLURM_CPUS_ON_NODE}

mkdir $OUTDIR $OUTDIR/failed

# trim r1 & r2 files
ls $INDIR/$FQPATTERN | \
	sed -e "s/$EXTPATTERN//" -e 's/.*\///g' | \
	uniq | \
	parallel --no-notice -j $THREADS \
	fastp \
		--in1 $INDIR/{}$FWDEXT \
		--in2 $INDIR/{}$REVEXT \
		--out1 $OUTDIR/{}fp2_r1.fq.gz \
		--out2 $OUTDIR/{}fp2_r2.fq.gz \
		--unpaired1 $OUTDIR/{}fp2_unprd.fq.gz \
		--unpaired2 $OUTDIR/{}fp2_unprd.fq.gz \
		--failed_out $OUTDIR/failed/{}fp2_fail.fq.gz \
		-h $OUTDIR/{}r1r2_fastp.html \
		-j $OUTDIR/{}r1r2_fastp.json \
		--detect_adapter_for_pe \
		--length_required $MINLEN \
		--cut_front \
		--cut_front_window_size 1 \
		--cut_front_mean_quality 20 \
		--cut_right \
		--cut_right_window_size 10 \
		--cut_right_mean_quality 20 \
		--disable_trim_poly_g \
		--correction \
		--disable_quality_filtering \
		--unqualified_percent_limit 40 \
		--report_title "Second Trim 4 De Novo R1R2"

# remove unprd reads

mv $OUTDIR/*unprd* $OUTDIR/failed

# run multiqc on the fp2 results

multiqc $OUTDIR -n $OUTDIR/2nd_fastp_report
