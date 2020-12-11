#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=25g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.filter_fastq_by_read_name
#PBS -j oe
#PBS -J 0-1
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1
MEMORY=20G

# input fastq file
IN_DIR=.
IN_SEQ=($(find ${IN_DIR} \( -name "*.fastq" -not -name "*converted.fastq" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.fastq}

# input read list
READ_LIST=$(find $IN_DIR -name "${BASE%.*}.converted_reads.txt")

# ----------------Commands------------------- #
# filter fastq by read names
filterbyname.sh -Xmx$MEMORY in=${FILE} out=${BASE}.converted.fastq names=${READ_LIST} include=t
