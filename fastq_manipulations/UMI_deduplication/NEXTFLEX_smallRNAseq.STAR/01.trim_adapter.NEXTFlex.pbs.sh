#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=30g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01a.trim_adapters.NF.bbduk
#PBS -J 0-1
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6

IN_DIR=../Links
IN_SEQ=($(find ${IN_DIR} -name "*.txt.gz"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.txt.gz}

ADAPTER=/common/WORK/fhorvat/Projekti/Svoboda/scripts/fastq_manipulations/adapter_trimming/NEXTflex_smallRNA_3prime_adapter.fa

BBDUK_PAR1="overwrite=t \
ktrim=r \
k=19 \
rcomp=t \
mink=10 \
hdist=1 \
minoverlap=8 \
threads=$THREADS"

BBDUK_PAR2="overwrite=t \
forcetrimright2=4 \
forcetrimleft=4 \
minlength=15 \
threads=$THREADS"

# ----------------Commands------------------- #
# trim right adapter
bbduk.sh -Xmx30g in=$FILE out=${BASE}.atrim.txt.gz ref=$ADAPTER stats=${BASE}.atrim.stats $BBDUK_PAR1

# force trim 4 bases from right and left
#bbduk.sh -Xmx30g in=${BASE}.atrim.txt.gz out=${BASE}.txt.gz stats=${BASE}.ftrim.stats $BBDUK_PAR2

# remove 
#[ -f "${BASE}.txt.gz" ] && rm ${BASE}.atrim.txt.gz
