#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.bam_read_list
#PBS -J 0-2
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
IN_DIR=.
IN_SEQ=($(find $IN_DIR -maxdepth 1 \( -name "*.bam" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.bam}

# ----------------Commands------------------- #
# get all unique reads from bam as a list
#sambamba view ${FILE} | awk '{ a[$1]++ } END { for (b in a) { print b } }' > ${BASE}.read_list.txt

# get all primary alignments read names from bam as a list
samtools view -F 256 ${FILE} | awk '{print $1}' > ${BASE}.read_list.txt
