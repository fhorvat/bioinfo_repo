#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02b.filter_bams.24to31
#PBS -l select=ncpus=1:mem=30g
#PBS -j oe
#PBS -J 0-3
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# source
source ./000.load_variables.sh

# files
INPUT_DIR=.
IN_SEQ=($(find $INPUT_DIR \( -name "*.bam" -not -name "*21to23nt*" -not -name "*24to31nt*" -not -name "*19to32*" -not -name "*perfect*" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}.24to31nt

# ----------------Commands------------------- #
# remove reads with deletions and insertions
# get reads between 24-31 nt long
samtools view -h ${FILE} | \
awk -F '\t' '($1 ~ /^@/ || $6 !~/I/) && ($1 ~ /^@/ || $6 !~/D/) && ($1 ~ /^@/ || $6 ~/24M|25M|26M|27M|28M|29M|30M|31M/)' | \
samtools view -Sb - > ${BASE}.bam

# index
samtools index ${BASE}.bam
