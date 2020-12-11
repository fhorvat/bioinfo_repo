#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=40g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.split_fastq_in_chunks
#PBS -j oe
#PBS -J 0-1
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1
MEMORY=100G

# input fastq file
IN_DIR=..
IN_SEQ=($(find ${IN_DIR} -maxdepth 1 \( -name "*.fastq" -not -name "*converted.fastq" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%fastq}

# ----------------Commands------------------- #
# split fastq into 10 million reads chunks
split -l 40000000 -d --additional-suffix=.fastq $FILE $BASE
