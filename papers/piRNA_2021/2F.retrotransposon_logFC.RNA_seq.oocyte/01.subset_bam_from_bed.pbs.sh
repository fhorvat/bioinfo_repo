#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.filter_bam
#PBS -l select=ncpus=6:mem=50g
#PBS -j oe
#PBS -J 0-6
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# number of threads
THREADS=6

# files
INPUT_DIR=/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.RNAseq/Data/Mapped/STAR_Siomi.multimappers
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 \( -name "*.bam" -not -name "*F66*" -not -name "*F53*" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}

BED=../../rmsk.Siomi.20200701.raw.SINE_LINE_LTRs.bed

# ----------------Commands------------------- #
# get only PERFECT reads overlapping .bed
#samtools view -@ ${THREADS} -b -L ${BED} ${FILE} > ./${BASE}.bam
samtools view -@ ${THREADS} -b -L ${BED} ${FILE} | \
bamtools filter -out ${BASE}.bam -tag "nM:0"

# index
samtools index -@ ${THREADS} ${BASE}.bam
