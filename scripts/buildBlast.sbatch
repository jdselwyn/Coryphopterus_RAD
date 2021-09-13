#!/bin/bash

#SBATCH --time=1-00:00:00
#SBATCH -p gpu,normal,cbirdq
#SBATCH --job-name=BlAST_Build
#SBATCH -o SLURM_out/buildBlast-%j.out

module load R/gcc/64/3.5.1
module load blast+/gcc7/2.9.0


# Step 1 - download all available mtDNA for any Coryphopterus
Rscript scripts/findCoryphopterusMito.R

#Step 2 - Create BLAST Database
#https://open.oregonstate.education/computationalbiology/chapter/command-line-blast/

makeblastdb -in Reference_Sequence/Coryphopterus_mtDNA.fasta \
  -out Reference_Sequence/CoryphopterusBlast \
  -dbtype nucl \
  -title "Coryphopterus Blast DB" \
  -parse_seqids
