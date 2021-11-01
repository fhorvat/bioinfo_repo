#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03.gzip_to_bgzip
#PBS -l select=ncpus=10:mem=5g
#PBS -J 0-8
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=10

INPUT_DIR=.
IN_SEQ=(`ls ${INPUT_DIR}/*.fa.gz`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fa.gz}

# ----------------Commands------------------- #
# ungzip, index, gzip
unpigz -p $THREADS ${FILE}
bgzip -@ $THREADS ${BASE}.fa
samtools faidx ${BASE}.fa.gz
