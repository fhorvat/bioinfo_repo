#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.trim_adapters.cutadapt
#PBS -J 0-6
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6

IN_DIR=.
IN_SEQ=($IN_DIR/s_GV_WT*.fq.gz)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.fq.gz}

# ----------------Commands------------------- #
# right trim
cutadapt -a CTGTCTCTTATA -a TATAAGAGACAG -q 20 --minimum-length=25 --cores=$THREADS -o ${BASE}.txt.gz $FILE
