#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N GarciaLopez_2015
#PBS -l select=ncpus=2:mem=10g
#PBS -J 0-4
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
REF_DIR=/common/WORK/fhorvat/reference/mouse/mm10/genome_indexes/STAR/sjdbOverhang_100

INPUT_DIR=..
IN_SEQ=($INPUT_DIR/*.bam)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE##*$INPUT_DIR/}
BASE=${BASE%_Aligned.sortedByCoord.out.bam}

MIN_LENGTH=21
MAX_LENGTH=23
BASE=${BASE}"_"$MIN_LENGTH"to"$MAX_LENGTH

# ----------------Commands------------------- #
# filter from .bam with alignment length from min to max
java -jar /common/WORK/fhorvat/programi/samjs/jvarkit/dist/samjs.jar -e '!record.readUnmappedFlag  && record.cigar.referenceLength  >= 21 && record.cigar.referenceLength  <= 23' $FILE | samtools view -Sb -o ${BASE}.bam -

# bam index
samtools index ${BASE}.bam

# bam to bedGraph, bedGraph to bigWig
genomeCoverageBed -ibam ${BASE}.bam -bg -split -g $REF_DIR/chrNameLength.txt > ${BASE}.bedGraph
wigToBigWig ${BASE}.bedGraph $REF_DIR/chrNameLength.txt ${BASE}.bw
[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph
