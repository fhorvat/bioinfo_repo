#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=30g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.trim_adapters.bbduk
#PBS -j oe
#PBS -J 0-5
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6

ADAPTER=../../Documentation/NEBNextSmallRNAadapter.fa
IN_DIR=../Links
IN_SEQ=($IN_DIR/*.txt.gz)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.txt.gz}

BBDUK_PAR="overwrite=t \
ktrim=r \
k=23 \
rcomp=t \
mink=12 \
hdist=1 \
minoverlap=8 \
minlength=15 \
threads=$THREADS"

# ----------------Commands------------------- #
bbduk.sh in=$FILE out=${BASE}_cl.txt.gz ref=$ADAPTER stats=${BASE}_cl.stats $BBDUK_PAR
