#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.scale_bigWig
#PBS -l select=ncpus=1:mem=10g
#PBS -J 0-5
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
STAR_INDEX=/common/DB/genome_reference/cow/bosTau9.ARS-UCD1.2.GCA_002263795.2/STAR_index.2.7/sjdbOverhang_100
CHR_LENGTH=${STAR_INDEX}/chrNameLength.txt

INPUT_DIR=.
IN_SEQ=(${INPUT_DIR}/*.bam)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}

LOG=../3_logs/log.read_stats.txt

# ----------------Commands------------------- #
# just in case
awk -v PAT="$BASE" -F '\t' '$1 ~ PAT {print $1}' $LOG

# get number of mapped reads in millions
RPM=`awk -v PAT="$BASE" -F '\t' '$1 ~ PAT {sum += $8} END {print sum}' $LOG`
RPM=`echo "scale=6; 1000000.0/"$RPM | bc`
BASE=${BASE}.scaled

# bam to scaled bedGraph, bedGraph to bigWig
genomeCoverageBed -ibam $FILE -bg -scale $RPM -split > ${BASE}.bedGraph
wigToBigWig ${BASE}.bedGraph $CHR_LENGTH ${BASE}.bw
[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph
