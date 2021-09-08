#!/bin/bash

#SBATCH --job-name=fltrVCF
#SBATCH -o SLURM_out/fltrVCF-%j.out
#SBATCH -p normal,cbirdq
#SBATCH --time=96:00:00
#SBATCH --nodes=1

#to run use the following command:
#sbatch scripts/fltrVCF.sbatch h_sp/fltVCF config.fltr.ind B

module load ddocent
module load perl/5.22.0
module load vcftools/0.1.15
module load rad_haplotyper/1.1.5
module load parallel
module load vcflib/1.0
module load samtools
module load R/gcc7/3.5.3

dir=$1
configFile=$2
prefixExtra=$3

mkdir -p ${dir}

dDocent="$(pwd)/scripts/dDocentHPC/dDocentHPC.bash"
fltrVCF="$(pwd)/scripts/fltrVCF/fltrVCF.bash"
rScripts="$(pwd)/scripts/fltrVCF/scripts"
HWE="$(pwd)/scripts/fltrVCF/filter_hwe_by_pop_HPC.pl"
haplotyper="$(pwd)/scripts/rad_haplotyper/rad_haplotyper.pl"

config="$(pwd)/${configFile}"
BAMS="$(pwd)/$(echo ${dir} | sed "s/fltVCF.*/mkVCF/g")"
BEDS="$(pwd)/$(echo ${dir} | sed "s/fltVCF.*/mkVCF/g")/mapped.2.2.bed"
startVCFS="$(pwd)/$(echo ${dir} | sed "s/fltVCF.*/mkVCF/g")/TotalRawSNPs.2.2.vcf"
reference="$(pwd)/$(echo ${dir} | sed "s/fltVCF.*/mkREF/g")/reference.2.2.fasta"
POPMAP="$(pwd)/$(echo ${dir} | sed "s/fltVCF.*/mkVCF/g")/popmap.2.2"
species="$(echo ${dir} | sed "s/\/fltVCF.*//g")_${prefixExtra}"

cd ${dir}

bash ${fltrVCF} -s ${config} \
  -b ${BAMS} \
  -d ${BEDS} \
  -v ${startVCFS} \
  -g ${reference} \
  -p ${POPMAP} \
  -o ${species} \
  -R ${rScripts} \
  -w ${HWE} \
  -r ${haplotyper} \
  -t ${SLURM_CPUS_ON_NODE}
