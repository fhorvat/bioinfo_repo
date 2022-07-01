#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.sort_bam
#PBS -l select=ncpus=12:mem=60g
#PBS -J 0-3
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# source
source ./000.load_variables.sh

# override number of threads 
THREADS=12

# input
INPUT_DIR=.
IN_SEQ=($(find ${INPUT_DIR} \( -name "*.bam" -not -name "*.sorted.bam" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}

# ----------------Commands------------------- #
# sort bam
samtools sort -@ ${THREADS} -o ${BASE}.sorted.bam ${BASE}.bam

# bam index
samtools index -@ ${THREADS} ${BASE}.sorted.bam
