#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.insert_size.bbmerge
#PBS -l select=ncpus=1:mem=5g
#PBS -J 0-29
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
INPUT_DIR=..
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 -name "*.txt.gz"))
UNIQ_SEQ=($(printf "%s\n" "${IN_SEQ[@]%_[1-2].txt.gz}" | sort -u))
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}

# ----------------Commands------------------- #
# get insert size with bbmerge
bbmerge.sh \
in1=${FILE}_1.txt.gz \
in2=${FILE}_2.txt.gz \
ihist=${BASE}.insertHis.txt \
outi=${BASE}.insertSize.txt \
loose \
reads=2m \
2> ${BASE}.log
