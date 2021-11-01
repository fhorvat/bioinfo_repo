#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.05.create_bigWig
#PBS -l select=ncpus=1:mem=10g
#PBS -J 0-1
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
CHR_LENGTH=./L1_1_CGr.1052918.chrNameLength.txt

IN_DIR=.
IN_SEQ=($(find ${IN_DIR} -maxdepth 1 \( -name "*.bam" -and -name "*merged*" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE%.merged.bam}

# ----------------Commands------------------- #
# bam to scaled bedGraph, bedGraph to bigWig
genomeCoverageBed -ibam $FILE -bg -split > ${BASE}.bedGraph
wigToBigWig ${BASE}.bedGraph $CHR_LENGTH ${BASE}.bw
[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph
