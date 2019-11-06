#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=30g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.trim_polyA.bbduk
#PBS -J 0-12
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6

IN_DIR=.
IN_SEQ=($(find $IN_DIR -name "*.txt.gz" -not -name "*polya*"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.txt.gz}

BBDUK_PAR="overwrite=t \
ktrim=r \
k=6 \
literal=AAAAAAAAA \
mm=f \
rcomp=t \
restrictright=40 \
hdist=1 \
threads=$THREADS"


# ----------------Commands------------------- #
# right trim
bbduk.sh in=$FILE out=${BASE}.polya.txt.gz stats=${BASE}.polya.stats $BBDUK_PAR
