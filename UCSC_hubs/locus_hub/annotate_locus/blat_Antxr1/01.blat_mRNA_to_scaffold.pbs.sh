#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.blat_to_lnc1
#PBS -l select=ncpus=1:mem=1g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1

INPUT_DIR=..
IN_SEQ=(${INPUT_DIR}/*.fa)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fa}

QUERY_INPUT_DIR=.
QUERY_FILE=(${QUERY_INPUT_DIR}/*.fa)
QUERY_BASE=${QUERY_FILE#${QUERY_INPUT_DIR}/}
QUERY_BASE=${QUERY_BASE%.fa}

# ----------------Commands------------------- #
# blat
blat -stepSize=5 -repMatch=1024 -minScore=0 -minIdentity=0 $FILE $QUERY_FILE ${BASE}.${QUERY_BASE}.all.psl

# get score 
pslScore ${BASE}.${QUERY_BASE}.all.psl > ${BASE}.${QUERY_BASE}.all.score.psl

# get .bed
#pslToBed ${BASE}.${QUERY_BASE}.psl ${BASE}.${QUERY_BASE}.bed
