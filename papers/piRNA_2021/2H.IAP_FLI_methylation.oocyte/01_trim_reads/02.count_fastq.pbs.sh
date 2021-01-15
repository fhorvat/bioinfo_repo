#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=60g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03.count_fastq
#PBS -j oe
#PBS -J 0-5
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1
MEMORY=60G

# input fastq file
IN_DIR=.
IN_SEQ=($(find ${IN_DIR} -maxdepth 1 \( -name "*.txt.gz" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.txt.gz}

# ----------------Commands------------------- #
# count fastq
COUNT=$(zcat $FILE | awk '{s++}END{print s/4}')

# print
echo -e ${BASE}"\t"${COUNT} >> fastq_files.counts.txt
