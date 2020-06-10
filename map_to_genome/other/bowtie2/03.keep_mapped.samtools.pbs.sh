#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=15g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.keep_mapped
#PBS -J 0-5
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6
MEMORY=15G

IN_DIR=.
IN_SEQ=($(find ${IN_DIR} -name "*bam" -not -name "*mapped*"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.bam}

# ----------------Commands------------------- #
# map with bowtie
samtools view -@$THREADS -Sb -F 4 $FILE > ${BASE}.mapped.bam
