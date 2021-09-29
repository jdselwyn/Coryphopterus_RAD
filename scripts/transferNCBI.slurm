#!/bin/bash

#SBATCH --job-name=uploadNCBI
#SBATCH -o SLURM_out/ncbi_upload_%j.out
#SBATCH -p normal,cbirdq
#SBATCH --time=96:00:00
#SBATCH --nodes=1

DIR=$1
NCBI_HOST=$2 #address from FTP upload option on: https://submit.ncbi.nlm.nih.gov/subs/sra/ (add s to front)
NCBI_USER=$3 #username from FTP upload option on: https://submit.ncbi.nlm.nih.gov/subs/sra/
NCBI_PW=$4 #password from FTP upload option on: https://submit.ncbi.nlm.nih.gov/subs/sra/
NCBI_FOLDER=$5 #folder they say to navigate to from FTP upload option on: https://submit.ncbi.nlm.nih.gov/subs/sra/
NCBI_FOLDER_MAKE=$6 #new folder for me to make - this is where everything will be

cd ${DIR}

sshpass -p ${NCBI_PW} sftp ${NCBI_USER}@${NCBI_HOST}

#to do after sftp success
cd ${NCBI_FOLDER}

mkdir ${NCBI_FOLDER_MAKE}
cd ${NCBI_FOLDER_MAKE}

mput *fq.gz
