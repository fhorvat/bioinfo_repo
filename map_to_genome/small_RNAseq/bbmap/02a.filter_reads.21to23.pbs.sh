#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02a.filter_bams.21to24
#PBS -l select=ncpus=1:mem=30g
#PBS -j oe
#PBS -J 0-9
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# names
STAR_INDEX=/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/STAR_index/sjdbOverhang_249
CHR_LENGTH=$STAR_INDEX/chrNameLength.txt

# files
INPUT_DIR=.
IN_SEQ=(`find $INPUT_DIR -name "*.bam" -not -name "*21to23nt*" -not -name "*24to31nt*" -not -name "*perfect*"`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}.21to23nt

# ----------------Commands------------------- #
# remove reads with deletions and insertions
# get reads between 21-23 nt long
# keep only perfect reads (nM tag == 0)
samtools view -h ${FILE} | \
awk -F '\t' '($1 ~ /^@/ || $6 !~/I/) && ($1 ~ /^@/ || $6 !~/D/) && ($1 ~ /^@/ || $6 ~/21=|22=|23=/)' | \
samtools view -Sb - > ${BASE}.bam

# index
samtools index ${BASE}.bam

# bam to bedGraph, bedGraph to bigWig
#genomeCoverageBed -ibam ${BASE}.bam -bg -split -g $CHR_LENGTH > ${BASE}.bedGraph
#wigToBigWig ${BASE}.bedGraph $CHR_LENGTH ${BASE}.bw
#[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph
