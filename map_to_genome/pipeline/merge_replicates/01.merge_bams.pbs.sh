#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.merge_bams
#PBS -l select=ncpus=6:mem=30g
#PBS -J 0-3
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------MAP MERGED------------------- #
# variables
THREADS=6

STAR_INDEX=/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/STAR_index/sjdbOverhang_249
CHR_LENGTH=${STAR_INDEX}/chrNameLength.txt

INPUT_DIR=..
IN_SEQ=(${INPUT_DIR}/*.bam)
STAGES=(${IN_SEQ[@]%_r[1-9]*.[S,P]E.bam})
STAGES=(${STAGES[@]%_old})
STAGES=(${STAGES[@]#${INPUT_DIR}/})
STAGES=($(printf "%s\n" ${STAGES[@]} | sort -u))

STAGE=${STAGES[$PBS_ARRAY_INDEX]}

FILES=($(find $INPUT_DIR -maxdepth 1 -name "$STAGE*r*.*E.bam"))
FILES=$(printf "%s " ${FILES[@]})

# just in case
echo $FILES

# ----------------Commands------------------- #
# merge
sambamba merge ${STAGE}.bam $FILES -t $THREADS

# bam to bedGraph, bedGraph to bigWig
genomeCoverageBed -ibam ${STAGE}.bam -bg -split -g $CHR_LENGTH > ${STAGE}.bedGraph
wigToBigWig ${STAGE}.bedGraph $CHR_LENGTH ${STAGE}.bw
[ -f "${STAGE}.bw" ] && rm ${STAGE}.bedGraph
