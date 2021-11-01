#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.scale_bigWig
#PBS -l select=ncpus=1:mem=20g
#PBS -J 0-%N_JOBS
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# source
source ./000.load_variables.sh

# input
INPUT_DIR=.
IN_SEQ=($(find ${INPUT_DIR} -name "*.bam"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}

# ----------------Commands------------------- #
# just in case
awk -v PAT="${BASE%${LAYOUT}}" -F '\t' '$1 ~ PAT {print $1}' ${LOG}

# get number of mapped reads in millions
SCALE_FACTOR=$(awk -v PAT="${BASE%${LAYOUT}}" -F '\t' '$1 ~ PAT {sum += $8} END {print 1000000/sum}' ${LOG})
BASE=${BASE}.scaled

# bam to scaled bedGraph, bedGraph to bigWig
genomeCoverageBed -ibam ${FILE} -bg -scale ${SCALE_FACTOR} -split > ${BASE}.bedGraph
wigToBigWig ${BASE}.bedGraph ${CHR_LENGTH} ${BASE}.bw
[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph
