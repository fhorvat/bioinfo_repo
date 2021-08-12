#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03.merge_bams
#PBS -l select=ncpus=10:mem=20g
#PBS -J 0-5
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------MAP MERGED------------------- #
# variables
THREADS=10
CHR_LENGTH=/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/STAR_index/sjdbOverhang_249/chrNameLength.txt

INPUT_DIR=.
IN_SEQ=(`find $INPUT_DIR -name "*genome.bam"`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE%.genome.bam}

# ----------------Commands------------------- #
# merge
sambamba merge -t $THREADS ${BASE}.bam ${BASE}.genome.bam ${BASE}.genome.merged.bam

# bam to bedGraph, bedGraph to bigWig
genomeCoverageBed -ibam ${BASE}.bam -bg -split -g $CHR_LENGTH > ${BASE}.bedGraph
wigToBigWig ${BASE}.bedGraph $CHR_LENGTH ${BASE}.bw
[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph
