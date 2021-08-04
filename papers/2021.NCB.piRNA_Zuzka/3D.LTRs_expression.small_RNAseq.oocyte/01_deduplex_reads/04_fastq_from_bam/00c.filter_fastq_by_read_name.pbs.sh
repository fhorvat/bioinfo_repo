#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.00c.filter_fastq_by_read_name
#PBS -j oe
#PBS -J 0-3
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1
MEMORY=20G

# input fastq file
IN_DIR=.
IN_SEQ=($(find ${IN_DIR} \( -name "*.fastq" -not -name "*filtered.fastq" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.fastq}

# input read list
READ_LIST=$(find $IN_DIR -name "${BASE}.18to32nt.read_list.txt")

# ----------------Commands------------------- #
# filter fastq by read names
filterbyname.sh -Xmx$MEMORY in=${FILE} out=${BASE}.filtered.fastq names=${READ_LIST} include=t
