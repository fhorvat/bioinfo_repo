#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.countFastq
#PBS -l select=ncpus=1:mem=20g
#PBS -J 0-15
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
INPUT_DIR=.
IN_SEQ=($INPUT_DIR/*txt.gz)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.txt.gz}

# ----------------Commands------------------- #
# counts reads in fastq
echo ${BASE} `zcat ${BASE}.txt.gz | awk '{s++}END{print s/4}'` >> stats.txt
# sort -t1 stats.txt > stats_sort.txt
