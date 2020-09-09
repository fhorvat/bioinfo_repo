#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=30g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.trim_adapters.RRBS.trim_galore
#PBS -J 0-5
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1

IN_DIR=../Links
IN_SEQ=($(find $IN_DIR -name "*.txt.gz" -not -name "*PE_*.txt.gz"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.txt.gz}

# ----------------Commands------------------- #
# trim
trim_galore --rrbs $FILE

# rename
mv ${BASE}.txt.gz_trimmed.fq.gz ${BASE}.txt.gz
