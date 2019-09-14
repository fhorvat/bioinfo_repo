#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.merge_fastq
#PBS -l select=ncpus=1:mem=30g
#PBS -J 0-4
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------MAP MERGED------------------- #
# variables
INPUT_DIR=.
IN_SEQ=($(find ${INPUT_DIR} -name "*.txt.gz" | grep -v "all"))
STAGES=($(printf "%s\n" "${IN_SEQ[@]%_r[1-9]*.[S,P]E.txt.gz}" | sort -u))
STAGES=(${STAGES[@]#${INPUT_DIR}/})
STAGE=${STAGES[$PBS_ARRAY_INDEX]}

FILES=($(find $INPUT_DIR -maxdepth 1 -name "$STAGE*r*.*E.txt.gz"))
FILES=$(printf "%s " ${FILES[@]})

# ----------------Commands------------------- #
# merge
cat $FILES > ${STAGE}_all.SE.txt.gz

