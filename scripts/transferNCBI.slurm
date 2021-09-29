#!/bin/bash

#SBATCH --job-name=uploadNCBI
#SBATCH -o SLURM_out/ncbi_upload_%j.out
#SBATCH -p normal,cbirdq
#SBATCH --time=96:00:00
#SBATCH --nodes=1

DIR=$1
NCBI_HOST=$2
NCBI_USER=$3
NCBI_PW=$4
NCBI_FOLDER=$5

cd ${DIR}

sshpass -p ${NCBI_PW} sftp ${NCBI_USER}@${NCBI_HOST}

#to do after sftp success
mkdir ${NCBI_FOLDER}
cd ${NCBI_FOLDER}

mput *fq.gz
