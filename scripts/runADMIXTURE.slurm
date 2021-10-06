#!/bin/bash

#SBATCH --job-name=admixture
#SBATCH -o SLURM_out/admixture-%j.out
#SBATCH -p normal,cbirdq
#SBATCH --time=96:00:00
#SBATCH --nodes=1

#https://speciationgenomics.github.io/ADMIXTURE/
#https://dalexander.github.io/admixture/admixture-manual.pdf

module load plink/1.90

dir=$1
inVCF=$2
maxK=$3
CV=$4

cutoffs=$(echo ${inVCF} | grep -Eo '\.[0-9]+\.[0-9]+\.' | grep -Eo '[0-9]+\.[0-9]+')
outNAME=$(echo ${inVCF} | grep -oP '(?<=/).*?(?=\.)')
outNAME="${outNAME}.${cutoffs}"

mkdir -p ${dir}

inVCF="$(pwd)/${inVCF}"
cd ${dir}

plink --vcf ${inVCF} --make-bed --maf 0.05 --recode --allow-extra-chr --out ${outNAME}
awk '{$1=0;print $0}' ${outNAME}.bim > ${outNAME}.bim.tmp
mv ${outNAME}.bim.tmp ${outNAME}.bim

for ((K=1; K<=${maxK}; K++))
do
 admixture --cv=${CV} $outNAME.bed $K -j${SLURM_CPUS_ON_NODE} > log${K}.out
done

grep -h CV log*.out
#Some R stuff here to make a ggplot of the CV


#Some R stuff here to make STRUCTURE plots for each value of K
