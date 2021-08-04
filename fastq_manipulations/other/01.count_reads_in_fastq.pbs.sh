#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.countFastq
#PBS -l select=ncpus=1:mem=20g
#PBS -J 0-3
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
INPUT_DIR=.
IN_SEQ=($(find ${INPUT_DIR} -name "*.fastq.gz" -or -name "*.txt.gz"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fastq.gz}
BASE=${BASE%.txt.gz}

# ----------------Commands------------------- #
# counts reads in fastq
echo ${BASE} `zcat ${FILE} | awk '{s++}END{print s/4}'` >> stats.txt
# sort -t1 stats.txt > stats_sort.txt
