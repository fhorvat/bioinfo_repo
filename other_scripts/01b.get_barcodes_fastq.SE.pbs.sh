#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.get_barcodes
#PBS -l select=ncpus=1:mem=10g
#PBS -J 0-8
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1
MEMORY=10g

INPUT_DIR=..
IN_SEQ=($(find $INPUT_DIR -name "*.txt.gz"))
UNIQ_SEQ=($(printf "%s\n" "${IN_SEQ[@]%.txt.gz}" | sort -u))
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}

# ----------------Commands------------------- #
zcat ${FILE}.txt.gz | awk '{if (NR % 4 == 1) print $2}' | awk -F ":" '{print $4}' | sort | uniq -c | sort -k1 -n -r > ${BASE}.barcodes.txt
