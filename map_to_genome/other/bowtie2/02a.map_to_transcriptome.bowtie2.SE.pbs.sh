#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=50g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.align.bowtie2
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

INDEX=bowtie_index

BOWTIE2_PAR="-p $THREADS \
-q \
--no-unal \
-x $INDEX"

# ----------------Commands------------------- #
# align reads to transcriptome
bowtie2 $BOWTIE2_PAR -1 ${FILE}_1.txt.gz -2 ${FILE}_2.txt.gz 2> ${BASE}.stats.txt | samtools view -@$THREADS -Sb -o ${BASE}.bam
