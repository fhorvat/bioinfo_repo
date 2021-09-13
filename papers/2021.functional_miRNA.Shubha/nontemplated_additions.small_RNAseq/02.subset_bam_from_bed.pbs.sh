#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.filter_bam
#PBS -l select=ncpus=1:mem=30g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# bed files
INPUT_DIR=..
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 -name "*.bed" -and -name "*pig*"))
FILE_BED=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE_BED#${INPUT_DIR}/}
BASE=${BASE%.top10_miRNA.bed}

# bam files
INPUT_DIR=../../datasets/${BASE}/5_perfect_reads/merged_replicates
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 \( -name "*.bam" -and \( -name "s_bosTau_GV.21to23nt.bam" -or -name "s_MII.21to23nt.bam" -or -name "s_oocyte_large.21to23nt.bam" \
-or -name "s_oocyte_small.21to23nt.bam" \) \)))
FILE_BAM=${IN_SEQ[1]}
BASE_BAM=${FILE_BAM#${INPUT_DIR}/}
BASE_BAM=${BASE_BAM%.bam}

# ----------------Commands------------------- #
# get reads overlapping bed
#bamtools filter -in $FILE -out ${BASE}
samtools view -b -L ${FILE_BED} ${FILE_BAM} > ${BASE}.top10_miRNA.${BASE_BAM}.bam

# index
samtools index ${BASE}.top10_miRNA.${BASE_BAM}.bam

# bam to bedGraph, bedGraph to bigWig
#genomeCoverageBed -ibam ${BASE}.bam -bg -split -g $CHR_LENGTH > ${BASE}.bedGraph
#wigToBigWig ${BASE}.bedGraph $CHR_LENGTH ${BASE}.bw
#[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph
