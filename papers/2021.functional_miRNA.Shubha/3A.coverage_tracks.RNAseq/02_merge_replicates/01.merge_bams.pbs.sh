#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.merge_bams
#PBS -l select=ncpus=6:mem=30g
#PBS -J 0-5
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------MAP MERGED------------------- #
# variables
THREADS=6

INPUT_DIR=..
IN_SEQ=(${INPUT_DIR}/*.bam)
STAGES=(${IN_SEQ[@]%.r[1-9]*.[S,P]E.*bam})
STAGES=(${STAGES[@]#${INPUT_DIR}/})
STAGES=($(printf "%s\n" ${STAGES[@]} | sort -u))

STAGE=${STAGES[$PBS_ARRAY_INDEX]}

FILES=($(find $INPUT_DIR -maxdepth 1 -name "$STAGE*r*.*E.*bam"))
FILES=$(printf "%s " ${FILES[@]})

# just in case
echo $FILES

# ----------------Commands------------------- #
# merge
sambamba merge ${STAGE}.bam $FILES -t $THREADS
