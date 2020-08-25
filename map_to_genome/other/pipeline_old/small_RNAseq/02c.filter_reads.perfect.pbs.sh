#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02c.filter_bams.perfect
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
BASE=${BASE%.bam}.perfect

# ----------------Commands------------------- #
# remove reads with deletions and insertions
# keep only perfect reads (nM tag == 0)
samtools view -h ${FILE} | \
awk -F '\t' '($1 ~ /^@/ || $6 !~/I/) && ($1 ~ /^@/ || $6 !~/D/)' | \
samtools view -Sb - | \
bamtools filter -out ${BASE}.bam -tag "nM:0"

# index
samtools index ${BASE}.bam

# bam to bedGraph, bedGraph to bigWig
#genomeCoverageBed -ibam ${BASE}.bam -bg -split -g $CHR_LENGTH > ${BASE}.bedGraph
#wigToBigWig ${BASE}.bedGraph $CHR_LENGTH ${BASE}.bw
#[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph
