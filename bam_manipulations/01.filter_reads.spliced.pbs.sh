#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.unspliced_reads
#PBS -l select=ncpus=1:mem=10g
#PBS -j oe
#PBS -J 0-1
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
INPUT_DIR=..
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 \( -name "*.bam" -not -name "*spliced.bam" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}

# ----------------Commands------------------- #
# filter spliced reads (CIGARs with N)
# for unspliced reads (CIGARs without Ns) use: '($1 ~ /^@/ || $6 !~ /N/)'
samtools view -h $FILE | \
awk -F '\t' '($1 ~ /^@/ || $6 ~ /N/)' | \
samtools view -Sb - > ${BASE}.spliced.bam
