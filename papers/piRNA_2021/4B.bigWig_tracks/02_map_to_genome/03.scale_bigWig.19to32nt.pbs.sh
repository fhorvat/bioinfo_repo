#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03.scale_bigWig
#PBS -l select=ncpus=1:mem=10g
#PBS -J 0-5
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
STAR_INDEX=/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed/STAR_index.2.7/sjdbOverhang_100
CHR_LENGTH=${STAR_INDEX}/chrNameLength.txt
LOG=./library_sizes.txt

INPUT_DIR=.
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 -name "*bam"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
PATTERN=${BASE%.*.bam}
BASE=${BASE%.bam}

# ----------------Commands------------------- #
# just in case
awk -v PAT="$PATTERN.*.19to32nt" -F '\t' '$1 ~ PAT {print $0}' $LOG

# get numbers of reads in millions
RPM=`awk -v PAT="$PATTERN.*.19to32nt" -F '\t' '$1 ~ PAT {sum += $2} END {print sum}' $LOG`
RPM=`echo "scale=6; 1000000.0/"$RPM | bc`
BASE=${BASE}.scaled

# bam to scaled bedGraph, bedGraph to bigWig
genomeCoverageBed -ibam $FILE -bg -scale $RPM -split > ${BASE}.bedGraph
wigToBigWig ${BASE}.bedGraph $CHR_LENGTH ${BASE}.bw
[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph
