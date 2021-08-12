#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.psl_to_bed
#PBS -l select=ncpus=1:mem=10g
#PBS -J 0-1
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1

INPUT_DIR=.
IN_SEQ=(`find ${INPUT_DIR} -name "*.psl" -not -name "*.score*"`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.psl}

# ----------------Commands------------------- #
# psl to bed
pslToBed $FILE ${BASE}.bed
