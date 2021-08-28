#!/bin/bash -l

#SBATCH --job-name=clmp_r12
#SBATCH -o SLURM_out/clmp_r1r2-%A_%a.out
#SBATCH -p normal,gpu,cbirdq
#SBATCH --time=24:00:00

module load bbmap/38.90
module load parallel

#export SINGULARITY_BIND=/home/e1garcia

#TEMPDIR=/scratch-lustre/cbird
#INDIR=/home/e1garcia/shotgun_PIRE/fq_4denovo_fp1
#OUTDIR=/home/e1garcia/shotgun_PIRE/fq_4denovo_fp1_clmp
FQPATTERN=${1}
TEMPDIR=/tmp
INDIR=${2}
OUTDIR=${3}
THREADS=1   #clumpify uses a ton of ram, be conservative
GROUPS=auto   #controls how much ram is used, refer to manual
RAMPERTHREAD=233g   #have had to set as high as 233g with groups=1

mkdir -p $OUTDIR

all_samples=$(ls $INDIR/$FQPATTERN | \
	sed -e 's/r1\.fq\.gz//' -e 's/.*\///g')
all_samples=($all_samples)

sample_name=${all_samples[${SLURM_ARRAY_TASK_ID}]}
echo ${sample_name}

clumpify.sh \
	in=$INDIR/${sample_name}r1.fq.gz \
	in2=$INDIR/${sample_name}r2.fq.gz \
	out=$OUTDIR/${sample_name}clmp.r1.fq.gz \
	out2=$OUTDIR/${sample_name}clmp.r2.fq.gz \
	groups=$GROUPS \
	overwrite=t \
	usetmpdir=t \
	tmpdir=$TEMPDIR \
	dedupe=t \
	addcount=t \
	subs=2 \
	containment=t \
	consensus=f \
	-Xmx$RAMPERTHREAD
