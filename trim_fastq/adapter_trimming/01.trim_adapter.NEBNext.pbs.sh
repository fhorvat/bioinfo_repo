#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=30g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.trim_adapters.NN.bbduk
#PBS -J 0-3
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6

IN_DIR=../Links
IN_SEQ=($IN_DIR/*.txt.gz)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.txt.gz}

ADAPTER=/common/WORK/fhorvat/Projekti/Svoboda/scripts/trim_fastq/adapter_trimming/NEBNext_smallRNA_3prime_adapter.fa
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
# right trim
bbduk.sh in=$FILE out=${BASE}.txt.gz ref=$ADAPTER stats=${BASE}.stats $BBDUK_PAR
