#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=2:mem=10g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N bbduk
#PBS -J 0-4
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
ADAPTER=/common/WORK/fhorvat/reference/adapters/Illumina_TruSeq_Small_RNA.fa 

IN_DIR=../Links
IN_SEQ=($IN_DIR/*.fastq.gz)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE##*$IN_DIR/}
BASE=${BASE%%.fastq.gz}

BBDUK_PAR="overwrite=t \
k=22 \
mink=12 \
hdist=1 \
threads=1 \
minlength=15 \
minoverlap=8"

# ----------------Commands------------------- #
bbduk2.sh in=$FILE rref=$ADAPTER out=${BASE}_cl.fastq.gz stats=${BASE}_cl.stats $BBDUK_PAR 2> ${BASE}_cl.log
