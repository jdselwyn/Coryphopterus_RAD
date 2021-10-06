#!/bin/bash

#SBATCH --job-name=admixture
#SBATCH -o SLURM_out/admixture-%j.out
#SBATCH -p normal,cbirdq
#SBATCH --time=96:00:00
#SBATCH --nodes=1

#https://speciationgenomics.github.io/ADMIXTURE/
#https://dalexander.github.io/admixture/admixture-manual.pdf

module load plink/1.90
module load admixture/1.3.0
module load R/gcc/64/3.5.1

dir=$1
inVCF=$2
maxK=$3
CV=$4

cutoffs=$(echo ${inVCF} | grep -Eo '\.[0-9]+\.[0-9]+\.' | grep -Eo '[0-9]+\.[0-9]+')
outNAME=$(echo ${inVCF} | grep -oP '(?<=/).*?(?=\.)')
outNAME="${outNAME}.${cutoffs}"

mkdir -p ${dir}

inVCF="$(pwd)/${inVCF}"
R_CV="$(pwd)/scripts/admixtureCVPlot.R"
R_STRUCTURE="$(pwd)/scripts/admixtureStructurePlot.R"
cd ${dir}

plink --vcf ${inVCF} --make-bed --recode --allow-extra-chr --out ${outNAME}
awk '{$1=0;print $0}' ${outNAME}.bim > ${outNAME}.bim.tmp
mv ${outNAME}.bim.tmp ${outNAME}.bim

for ((K=1; K<=${maxK}; K++))
do
 admixture --cv=${CV} $outNAME.bed $K -j${SLURM_CPUS_ON_NODE} > ${outNAME}.log${K}.out
done

grep -h CV ${outNAME}.log*.out > ${outNAME}.CV_error
minK=$(Rscript ${R_CV} ${outNAME})
echo "Best K = ${minK}"

Rscript ${R_STRUCTURE} ${outNAME} ${minK}
