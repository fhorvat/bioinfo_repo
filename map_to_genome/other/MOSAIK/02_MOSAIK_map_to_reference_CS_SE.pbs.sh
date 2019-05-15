#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N MOSAIK_CS
#PBS -J 0-1
#PBS -l select=ncpus=6:mem=50g
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=4

REF_DIR=/common/WORK/fhorvat/reference/pig/susScr3/genome_indexes/MOSAIK
CHROM_SIZES=/common/WORK/fhorvat/reference/pig/susScr3/UCSC/susScr3.chrom.sizes

INPUT_DIR=/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/Cao_2014_pig/Data/Mapped/MOSAIK
IN_SEQ=($INPUT_DIR/*_cs.dat)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE##*$INPUT_DIR/}
BASE=${BASE%_cs.dat}

MOSAIK_PAR="-ibs ${REF_DIR}/mosaik_susScr3_bs.dat \
-ia ${REF_DIR}/mosaik_susScr3_cs.dat \
-j ${REF_DIR}/mosaik_susScr3_10 \
-annse /common/WORK/fhorvat/programi/MOSAIK/src/networkFile/2.1.78.se.ann \
-annpe /common/WORK/fhorvat/programi/MOSAIK/src/networkFile/2.1.78.pe.ann \
-hs 10 \
-mm 4 \
-mhp 100 \
-act 25 \
-p $THREADS"

# ----------------Commands------------------- #
# Map SOLiD reads to consensus
MosaikAligner -in $FILE -out ${BASE}_aligned.dat $MOSAIK_PAR &> ${BASE}_aligned.log

# convert to .bam
MosaikText -in ${BASE}_aligned.dat -bam ${BASE}
#[ -f "${BASE}.bam" ] && rm -f ${BASE}_aligned.dat

# sort, index
samtools sort -@ $THREADS ${BASE}_aligned.dat.bam ${BASE}_sorted
samtools index ${BASE}_sorted.bam
#[ -f "${BASE}_sorted.bam" ] && rm -f ${BASE}.bam

# bigwig
genomeCoverageBed -ibam ${BASE}_sorted.bam -bg -split -g $CHROM_SIZES > ${BASE}.bedGraph
wigToBigWig ${BASE}.bedGraph $CHROM_SIZES ${BASE}.bw
#[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph
