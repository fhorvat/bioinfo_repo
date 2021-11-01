#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=10:mem=60g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.trim_adapters.bbduk
#PBS -J 0-22
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=10

IN_DIR=../Links
IN_SEQ=(`ls ${IN_DIR}/*.txt.gz | grep "bosTau_"`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.txt.gz}

ADAPTER1=../../Documentation/adapter.fa
BBDUK_PAR1="overwrite=t \
ktrim=r \
k=8 \
rcomp=t \
mink=6 \
hdist=1 \
minoverlap=6 \
minlength=15 \
threads=$THREADS"

BBDUK_PAR2="overwrite=t \
forcetrimright2=4 \
forcetrimleft=4 \
minlength=15 \
threads=$THREADS"

# ----------------Commands------------------- #
# right trim
bbduk.sh in=$FILE out=${BASE}.trim1.txt.gz ref=$ADAPTER1 stats=${BASE}.trim1.stats $BBDUK_PAR1

# trim 4N bases from right and left
bbduk.sh in=${BASE}.trim1.txt.gz out=${BASE}.trim2.txt.gz stats=${BASE}.trim2.stats $BBDUK_PAR2

