#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.merge_with_extending
#PBS -l select=ncpus=10:mem=125g
#PBS -J 0-2
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=10
MEMORY=120g

ADAPTERS=../../../Documentation/bbmerge_adapters.fasta

INPUT_DIR=../../Links
IN_SEQ=(`ls ${INPUT_DIR}/*.txt.gz`)
UNIQ_SEQ=(`printf "%s\n" "${IN_SEQ[@]%_*txt.gz}" | sort -u`)
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}.merged

BBMERGE_PAR="adapter=$ADAPTERS \
rem \
k=62 \
extend2=50 \
ecct \
vstrict \
t=$THREADS \
prealloc=t \
-Xmx$MEMORY"

# ----------------Commands------------------- #
# merging
bbmerge.sh in1=${FILE}_1.txt.gz in2=${FILE}_2.txt.gz out=${BASE}.txt.gz $BBMERGE_PAR
