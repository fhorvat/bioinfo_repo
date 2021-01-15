#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.06.scale_bigWig
#PBS -l select=ncpus=1:mem=10g
#PBS -J 0-8
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
CHR_LENGTH=/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1/STAR_index/sjdbOverhang_100/chrNameLength.txt

INPUT_DIR=.
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 -name "*.bam" -not -name "*genome*" -not -name "*rDNA_45S*" -not -name "*sortedByName*"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}

# ----------------Commands------------------- #
# get number of mapped reads in millions
RPM=`grep ${BASE} log.read_stats.txt | awk -F "\t" '{print $8}'`
RPM=`echo "scale=6; 1000000.0/"$RPM | bc`
BASE=${BASE}.scaled

# bam to scaled bedGraph, bedGraph to bigWig
genomeCoverageBed -ibam $FILE -bg -scale $RPM -split -g $CHR_LENGTH > ${BASE}.bedGraph
wigToBigWig ${BASE}.bedGraph $CHR_LENGTH ${BASE}.bw
[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph
