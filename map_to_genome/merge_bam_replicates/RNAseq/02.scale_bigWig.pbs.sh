#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.scale_bigWig
#PBS -l select=ncpus=1:mem=40g
#PBS -J 0-%N_REPLICATES
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

LOG=../3_logs/log.read_stats.txt

# ----------------Commands------------------- #
# just in case
awk -v PAT="$BASE" -F '\t' '$1 ~ PAT {print $1}' $LOG

# get number of mapped reads in millions
SCALE_FACTOR=$(awk -v PAT="$BASE" -F '\t' '$1 ~ PAT {sum += $8} END {print 1000000/sum}' $LOG)
BASE=${BASE}.scaled

# bam to scaled bedGraph, bedGraph to bigWig
genomeCoverageBed -ibam $FILE -bg -scale $RPM -split -g $CHR_LENGTH > ${BASE}.bedGraph
wigToBigWig ${BASE}.bedGraph $CHR_LENGTH ${BASE}.bw
[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph
