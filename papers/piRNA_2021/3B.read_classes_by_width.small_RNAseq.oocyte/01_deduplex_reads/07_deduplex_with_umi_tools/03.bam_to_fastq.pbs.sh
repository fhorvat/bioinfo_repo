#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03.bam_to_fastq
#PBS -J 0-2
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12

IN_DIR=.
IN_SEQ=($(find $IN_DIR -maxdepth 1 \( -name "*.bam" -and -name "*18to32nt.bam" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.18to32nt.bam}

# ----------------Commands------------------- #
# get mapped reads to fastq
#samtools fastq -@ $THREADS ${FILE} > ${BASE}.fastq

# gzip
pigz -p $THREADS ${BASE}.fastq

# rename
mv ${BASE}.fastq.gz ${BASE}.txt.gz
