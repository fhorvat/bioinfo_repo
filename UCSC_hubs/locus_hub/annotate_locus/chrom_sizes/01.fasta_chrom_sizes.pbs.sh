#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.fasta_chrom_sizes
#PBS -l select=ncpus=1:mem=1g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1

INPUT_DIR=..
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 -name "*.fa.fai"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fa.fai}

# ----------------Commands------------------- #
# fasta index to chromosome sizes
cut -f1,2 $FILE > ${BASE}.chrom.sizes
