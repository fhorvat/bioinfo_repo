#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.blat_ovomucin
#PBS -l select=ncpus=1:mem=10g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1

INPUT_DIR=../..
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 -name "*.fa"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fa}

QUERY_INPUT_DIR=../flanking_genes.mm10
QUERY_FILE=${QUERY_INPUT_DIR}/Gfpt1.mRNA.fa
QUERY_BASE=${QUERY_FILE#${QUERY_INPUT_DIR}/}
QUERY_BASE=${QUERY_BASE%.fa}

# ----------------Commands------------------- #
# blat
blat -stepSize=5 -repMatch=1024 -minScore=0 -minIdentity=0 $FILE $QUERY_FILE ${BASE}.${QUERY_BASE}.psl

# get score 
pslScore ${BASE}.${QUERY_BASE}.psl > ${BASE}.${QUERY_BASE}.score.psl
