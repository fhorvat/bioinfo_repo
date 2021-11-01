#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.fasta_to_2bit
#PBS -l select=ncpus=1:mem=1g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1

INPUT_DIR=.
IN_SEQ=(`find ${INPUT_DIR} -name "*.fa"`)
FILE=${IN_SEQ[0]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fa}

# ----------------Commands------------------- #
# index fasta file
samtools faidx $FILE

# fasta to 2bit
faToTwoBit $FILE ${BASE}.2bit
