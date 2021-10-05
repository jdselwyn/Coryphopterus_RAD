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
module load R/gcc/64/3.5.1

dir=$1
startVCFS=$2
configFile=$3
prefixExtra=$4

mkdir -p ${dir}

dDocent="$(pwd)/scripts/dDocentHPC/dDocentHPC.bash"
fltrVCF="$(pwd)/scripts/fltrVCF/fltrVCF.bash"
rScripts="$(pwd)/scripts/fltrVCF/scripts"
HWE="$(pwd)/scripts/fltrVCF/filter_hwe_by_pop_HPC.pl"
haplotyper="$(pwd)/scripts/rad_haplotyper/rad_haplotyper.pl"


cutoffs=$(echo ${startVCFS} | grep -Eo '\.[0-9]+\.[0-9]+\.' | grep -Eo '[0-9]+\.[0-9]+')
helperVCF=$(echo ${startVCFS} | grep -oP '(?<=_).*?(?=/)')


config="$(pwd)/${configFile}"
BAMS="$(pwd)/$(echo ${dir} | sed "s/fltrVCF.*/mkVCF_${helperVCF}/g")"
BEDS="$(pwd)/$(echo ${dir} | sed "s/fltrVCF.*/mkVCF_${helperVCF}/g")/mapped.${cutoffs}.bed"
# startVCFS="$(pwd)/$(echo ${dir} | sed "s/fltVCF.*/mkVCF/g")/TotalRawSNPs.2.2.vcf"
startVCFS="$(pwd)/${startVCFS}"

if [ ${helperVCF} = "test" ]; then
	reference="$(pwd)/$(echo ${dir} | sed "s/fltrVCF.*/mkREF_MiSeq/g")/reference.${cutoffs}.fasta" #if test make MiSeq
else
	reference="$(pwd)/$(echo ${dir} | sed "s/fltrVCF.*/mkREF_${helperVCF}/g")/reference.${cutoffs}.fasta"
fi

POPMAP="$(pwd)/$(echo ${dir} | sed "s/fltrVCF.*/mkVCF_${helperVCF}/g")/popmap.${cutoffs}"
species="$(echo ${dir} | sed "s/fltrVCF.*/${helperVCF}/g")_${prefixExtra}"

cd ${dir}

bash ${fltrVCF} -s ${config} \
  -c ${cutoffs} \
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


#Summary Plot
bash ${rScripts}/fltrVCFstats2.sbatch ${SLURM_CPUS_ON_NODE} ${species} ${rScripts}
