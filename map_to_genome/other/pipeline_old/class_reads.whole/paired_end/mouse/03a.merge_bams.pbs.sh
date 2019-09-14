#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03a.merge_bams
#PBS -l select=ncpus=1:mem=60g
#PBS -J 0-%N_SAMPLES
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------MAP MERGED------------------- #
# variables
CHR_LENGTH=%GENOME_DIR/STAR_index/%SJDB_OVERHANG/chrNameLength.txt

INPUT_DIR=.
IN_SEQ=(`ls *.genome.Aligned.sortedByCoord.out.bam | grep -v "merged"`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE%.genome.Aligned.sortedByCoord.out.bam}.total

# ----------------Commands------------------- #
# merge
samtools merge ${BASE}.bam ${BASE%.total}.genome.Aligned.sortedByCoord.out.bam ${BASE%.total}.genome.merged.Aligned.sortedByCoord.out.bam

# index
samtools index ${BASE}.bam

# bam to bedGraph, bedGraph to bigWig
genomeCoverageBed -ibam ${BASE}.bam -bg -split -g $CHR_LENGTH > ${BASE}.bedGraph
wigToBigWig ${BASE}.bedGraph $CHR_LENGTH ${BASE}.bw
[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph
