#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.06.scale_bigWig
#PBS -l select=ncpus=1:mem=10g
#PBS -J 0-%N_SAMPLES
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
CHR_LENGTH=%GENOME_DIR/STAR_index/%SJDB_OVERHANG/chrNameLength.txt

INPUT_DIR=.
IN_BAM=(`ls ${INPUT_DIR}/*.genome.Aligned.sortedByCoord.out.bam`)
FILE=${IN_BAM[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.Aligned.sortedByCoord.out.bam}

# ----------------Commands------------------- #
# get number of mapped reads in millions
RPM=`grep ${BASE%.genome} log.read_stats.txt | awk -F "\t" '{print $8}'`
RPM=`echo "scale=6; 1000000.0/"$RPM | bc`
BASE=${BASE}".scaled.RPM_minus_rDNA"

# bam to scaled bedGraph, bedGraph to bigWig
genomeCoverageBed -ibam $FILE -bg -scale $RPM -split -g $CHR_LENGTH > ${BASE}.bedGraph
wigToBigWig ${BASE}.bedGraph $CHR_LENGTH ${BASE}.bw
[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph
