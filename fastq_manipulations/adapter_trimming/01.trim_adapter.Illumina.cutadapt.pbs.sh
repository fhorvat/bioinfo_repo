#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=30g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.trim_adapters.Illumina.cutadapt
#PBS -J 0-1
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6

IN_DIR=../Links
IN_SEQ=($IN_DIR/*.txt.gz)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.txt.gz}

# ----------------Commands------------------- #
# trim
cutadapt -a TGGAATTCTCGGGTGCCAAGG -m 18 -j $THREADS -o ${BASE}.txt.gz $FILE
