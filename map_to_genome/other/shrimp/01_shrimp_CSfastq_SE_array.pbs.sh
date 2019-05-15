#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N shrimpCS
#PBS -J 0-7
#PBS -l select=ncpus=6:mem=80g
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=4
REF_DIR=/common/WORK/fhorvat/reference/pig/susScr3
REF_GENOME=${REF_DIR}/genome_indexes/shrimp/susScr3-cs
REF_CHR_SIZE=${REF_DIR}/UCSC/susScr3.chrom.sizes

INPUT_DIR=/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/Cao_2014_pig/Data/Raw/Fastq
IN_SEQ=($INPUT_DIR/*.fastq.gz)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE##*$INPUT_DIR/}
BASE=${BASE%.fastq.gz}

SHRIMP_PAR="--threads $THREADS \
-L $REF_GENOME \
-o 1 \
-h 80% \
-E"

# ----------------Commands------------------- #
# mapping
gmapper-cs ${FILE} $SHRIMP_PAR > ${BASE}.sam 2> ${BASE}.log

# sam to bam
samtools view ${BASE}.sam -@ $THREADS -Sb > ${BASE}.bam
samtools sort -@ $THREADS ${BASE}.bam ${BASE}_sorted
samtools index ${BASE}_sorted.bam
[ -f "${BASE}_sorted.bam" ] && rm -f ${BASE}.bam
[ -f "${BASE}_sorted.bam" ] && rm -f ${BASE}.sam

# bam to bedGraph,  bedGraph to bigWig
genomeCoverageBed -ibam ${BASE}_sorted.bam -bg -split -g $REF_CHR_SIZE > ${BASE}.bedGraph
wigToBigWig ${BASE}.bedGraph $REF_CHR_SIZE ${BASE}.bw
[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph
