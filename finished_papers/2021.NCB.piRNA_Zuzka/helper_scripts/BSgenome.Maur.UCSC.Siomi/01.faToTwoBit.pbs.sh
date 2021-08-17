#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.faToTwoBit
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script parameters
THREADS=1

# set input variables
IN_DIR=.
IN_GENOME=($(find $IN_DIR -name "*.fasta"))
BASE=${IN_GENOME#${IN_DIR}/}
BASE=${BASE%.fasta}

# ----------------Commands------------------- #
# convert
faToTwoBit ${IN_GENOME} ${BASE}.2bit
