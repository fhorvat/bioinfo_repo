#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=50g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.align.salmon
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=50G

BASE_DIR=`basename ${PWD/trinity*/}`
INPUT_DIR=../../../../Raw/$BASE_DIR/Cleaned
IN_SEQ=(`ls ${INPUT_DIR}/s_*all.PE_[1,2].txt.gz`)
UNIQ_SEQ=(`printf "%s\n" "${IN_SEQ[@]%_*txt.gz}" | sort -u`)
FILE=${UNIQ_SEQ[0]}
BASE=${FILE#${INPUT_DIR}/}

INDEX=./salmon_index

# ----------------Commands------------------- #
# align reads to transcriptome
salmon quant -i $INDEX -l A -1 ${FILE}_1.txt.gz -2 ${FILE}_2.txt.gz --validateMappings -o transcripts_quant -p $THREADS 2> ${BASE}.stats.txt 
