#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.illumina_download
#PBS -l select=ncpus=1:mem=20g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# token generated from website here: /home/students/fhorvat/.basespace
# change username with: bs auth --force

# set project name
PROJECT_NAME="20190918TaborskaGVMILI"

# ----------------Commands------------------- #
# download data
bs download project -v -n $PROJECT_NAME -o $PROJECT_NAME
