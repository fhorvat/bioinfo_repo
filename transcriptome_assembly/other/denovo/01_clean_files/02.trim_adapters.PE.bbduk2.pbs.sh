#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=30g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.trim_adapters.PE
#PBS -j oe
#PBS -J 0-2
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6
MEMORY=30g

IN_DIR=..
IN_SEQ=($IN_DIR/s_GV*.txt.gz)
UNIQ_SEQ=(`printf "%s\n" "${IN_SEQ[@]%_*txt.gz}" | sort -u`)
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}

ADAPTER=./s_GV_WT_r1.PE.adapters.fasta
BBDUK_PAR="overwrite=t \
ktrim=r \
k=23 \
rcomp=t \
mink=12 \
hdist=1 \
minoverlap=8 \
minlength=25 \
tbo \
threads=$THREADS \
-Xmx$MEMORY"

# ----------------Commands------------------- #
# right trim
bbduk.sh in1=${FILE}_1.txt.gz in2=${FILE}_2.txt.gz out1=${BASE}_1.trim.txt.gz out2=${BASE}_2.trim.txt.gz ref=$ADAPTER stats=${BASE}.stats $BBDUK_PAR 2> ${BASE}.trim.log
