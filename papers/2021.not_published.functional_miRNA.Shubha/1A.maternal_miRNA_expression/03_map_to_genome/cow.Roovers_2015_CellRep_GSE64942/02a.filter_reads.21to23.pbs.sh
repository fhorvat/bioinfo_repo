#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02a.filter_bams.21to23
#PBS -l select=ncpus=1:mem=30g
#PBS -j oe
#PBS -J 0-5
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# source
source ./000.load_variables.sh

# files
INPUT_DIR=.
IN_SEQ=($(find $INPUT_DIR \( -name "*.bam" -not -name "*21to23nt*" -not -name "*24to31nt*" -not -name "*19to32*" -not -name "*perfect*" -and -name "*MI*" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}.21to23nt

# ----------------Commands------------------- #
# remove reads with deletions and insertions
# get reads between 21-23 nt long
samtools view -h ${FILE} | \
awk -F '\t' '($1 ~ /^@/ || $6 !~/I/) && ($1 ~ /^@/ || $6 !~/D/) && ($1 ~ /^@/ || $6 ~/21M|22M|23M/)' | \
samtools view -Sb - > ${BASE}.bam

# index
samtools index ${BASE}.bam
