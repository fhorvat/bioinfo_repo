#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03.filter_bams
#PBS -l select=ncpus=1:mem=30g
#PBS -j oe
#PBS -J 0-15
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# names
STAR_INDEX=/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/STAR_index/sjdbOverhang_249
CHR_LENGTH=$STAR_INDEX/chrNameLength.txt

# files
INPUT_DIR=.
IN_SEQ=(`find $INPUT_DIR -name "*.bam" | grep -v "21to23nt"`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}.21to23nt

# ----------------Commands------------------- #
# remove reads with deletions and insertions
# get reads between 21-23 nt long
# remove reads with more than one soft-clipped base at the 5' end
# get only perfect reads (nM tag == 0)
samtools view -h ${FILE} | \
awk -F '\t' '($1 ~ /^@/ || $6 !~/I/) && ($1 ~ /^@/ || $6 !~/D/) && ($1 ~ /^@/ || $6 ~/21M|22M|23M/)' | \
awk '{S=0; split($6,C,/[0-9]*/); n=split($6,L,/[NMSID]/);  if (and($2,0x10)>0 && C[n]=="S") {S=L[n-1]} else if (and($2,0x10)==0 && C[2]=="S") {S=L[1]}; if (S<=1) print }' | \
samtools view -Sb - | \
bamtools filter -out ${BASE}.bam -tag "nM:0"

# index
samtools index ${BASE}.bam

# bam to bedGraph, bedGraph to bigWig
genomeCoverageBed -ibam ${BASE}.bam -bg -split -g $CHR_LENGTH > ${BASE}.bedGraph
wigToBigWig ${BASE}.bedGraph $CHR_LENGTH ${BASE}.bw
[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph
