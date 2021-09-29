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
pattern='*fq.gz'
cd ${DIR}

#make file with commands of what to do after sftp connection
echo cd ${NCBI_FOLDER} > sftp_commands.txt
echo mkdir ${NCBI_FOLDER_MAKE} >> sftp_commands.txt
echo cd ${NCBI_FOLDER_MAKE} >> sftp_commands.txt
printf 'mput *fq.gz\n' >> sftp_commands.txt

sshpass -p ${NCBI_PW} sftp -b sftp_commands.txt ${NCBI_USER}@${NCBI_HOST}

rm sftp_commands.txt
