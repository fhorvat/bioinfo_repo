#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.scaffold_from_fasta
#PBS -l select=ncpus=1:mem=1g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1

INPUT_DIR=..
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 -name "*.fa"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fa}

SCAFFOLD="PVKX01000823.1"

# ----------------Commands------------------- #
# get fasta on correct strand with bedtools
samtools faidx ${FILE} ${SCAFFOLD} > ${BASE}.${SCAFFOLD}.fa
