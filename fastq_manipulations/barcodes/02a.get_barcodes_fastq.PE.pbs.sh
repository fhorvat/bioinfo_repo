#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.get_barcodes
#PBS -l select=ncpus=1:mem=20g
#PBS -J 0-11
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1
MEMORY=20g

INPUT_DIR=..
IN_SEQ=($(find $INPUT_DIR -name "*.txt.gz"))
UNIQ_SEQ=($(printf "%s\n" "${IN_SEQ[@]%_*txt.gz}" | sort -u))
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}

# ----------------Commands------------------- #
zcat ${FILE}_1.txt.gz | awk '{if (NR % 4 == 1) print $2}' | awk -F ":" '{print $4}' | sort | uniq -c | sort -k1 -n -r > ${BASE}_1.barcodes.txt
zcat ${FILE}_2.txt.gz | awk '{if (NR % 4 == 1) print $2}' | awk -F ":" '{print $4}' | sort | uniq -c | sort -k1 -n -r > ${BASE}_2.barcodes.txt
