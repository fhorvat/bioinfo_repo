#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.get_adapters
#PBS -l select=ncpus=6:mem=40g
#PBS -J 0-2
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6
MEMORY=40g

INPUT_DIR=..
IN_SEQ=(`ls ${INPUT_DIR}/s_GV*.txt.gz`)
UNIQ_SEQ=(`printf "%s\n" "${IN_SEQ[@]%_*txt.gz}" | sort -u`)
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}

# ----------------Commands------------------- #
# merging
bbmerge.sh in1=${FILE}_1.txt.gz in2=${FILE}_2.txt.gz outa=${BASE}.adapters.fasta t=$THREADS -Xmx$MEMORY
