#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=30g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.fastq_to_fasta
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
DIR_BASE=`basename ${PWD/oases/}`

INPUT_DIR=../../../Raw/$DIR_BASE/Cleaned
IN_SEQ=(`ls ${INPUT_DIR}/*all.PE_[1,2].txt.gz`)
UNIQ_SEQ=(`printf "%s\n" "${IN_SEQ[@]%_[1-2].txt.gz}" | sort -u`)
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.txt.gz}

# ----------------Commands------------------- #
# get insert size
bbmerge.sh in1=${FILE}_1.txt.gz in2=${FILE}_2.txt.gz ihist=${BASE}.hist.txt reads=2m
