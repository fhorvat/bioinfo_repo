#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=30g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.trim_adapters.Truseq.bbduk
#PBS -J 0-10
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6
MEMORY=30g

IN_DIR=../Links
IN_SEQ=($(find $IN_DIR -name "*.txt.gz"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.txt.gz}

ADAPTER=/common/WORK/fhorvat/Projekti/Svoboda/scripts/trim_fastq/adapter_trimming/Illumina_TruSeq_Small_RNA.fa

BBDUK_PAR="overwrite=t \
ktrim=r \
k=22 \
rcomp=t \
mink=12 \
hdist=1 \
minoverlap=8 \
minlength=15 \
threads=$THREADS"

# ----------------Commands------------------- #
# right trim
bbduk.sh -Xmx$MEMORY in=$FILE out=${BASE}.txt.gz ref=$ADAPTER stats=${BASE}.stats $BBDUK_PAR
