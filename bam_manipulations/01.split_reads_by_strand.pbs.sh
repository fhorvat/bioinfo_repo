#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.split_bam
#PBS -l select=ncpus=1:mem=30g
#PBS -J 0-11
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------MAP MERGED------------------- #
# variables
THREADS=1

STAR_INDEX=/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/STAR_index/sjdbOverhang_249
CHR_LENGTH=${STAR_INDEX}/chrNameLength.txt

INPUT_DIR=..
FILES=($(find $INPUT_DIR -maxdepth 1 -name "*.bam"))
FILE=${FILES[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}

# ----------------Commands------------------- #
# bamtools filter by strand
bamtools filter -in $FILE -out ${BASE}.minus.bam -isReverseStrand true
bamtools filter -in $FILE -out ${BASE}.plus.bam -isReverseStrand false
