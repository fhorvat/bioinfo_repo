#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=10g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.04.remove_chromosome
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1
MEMORY=10G

IN_DIR=.
IN_SEQ=($(find ${IN_DIR} -maxdepth 1 -name "*merged.bam"))
FILE=${IN_SEQ[0]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.bam}

# ----------------Commands------------------- #
# remove chromosme with "/" from bam file
samtools idxstats $FILE | cut -f 1 | grep -v "/" | xargs samtools view -b $FILE > ${BASE}.filt.bam
