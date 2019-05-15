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
INPUT_DIR=../../../../Raw/$BASE_DIR/Cleaned/extended_merge.bbnorm_normalized
IN_SEQ=(`ls ${INPUT_DIR}/s_*all.PE.txt.gz`)
FILE=${IN_SEQ[0]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.PE.txt.gz}

INDEX=bowtie_index

BOWTIE2_PAR="-p $THREADS \
-q \
--no-unal \
-x $INDEX"

# ----------------Commands------------------- #
# align reads to transcriptome
bowtie2 $BOWTIE2_PAR -r ${FILE} 2> ${BASE}.SE.stats.txt | samtools view -@$THREADS -Sb -o ${BASE}.SE.bam
