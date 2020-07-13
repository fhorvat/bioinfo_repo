#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=30g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.merge_mapped
#PBS -J 0-1
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6
MEMORY=30G

PATTERNS=(s_mesAur_3kb s_mesAur_fragment)
PATTERN=${PATTERNS[$PBS_ARRAY_INDEX]}

IN_DIR=.
IN_SEQ=($(find ${IN_DIR} -maxdepth 1 \( -name "${PATTERN}*" -and -name "*.bam" -not -name "*merged*" \)))
FILES=$(printf "%s " ${IN_SEQ[@]})
BASE=${PATTERN}

# ----------------Commands------------------- #
# merge all mapped files
samtools merge -@$THREADS ${BASE}.merged.bam ${FILES}

# index
samtools index ${BASE}.merged.bam
