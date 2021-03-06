#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.05.scale_bigWig
#PBS -l select=ncpus=1:mem=10g
#PBS -J 0-31
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
STAR_INDEX=/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/STAR_index/sjdbOverhang_249
CHR_LENGTH=$STAR_INDEX/chrNameLength.txt

INPUT_DIR=.
IN_BAM=(${INPUT_DIR}/*.bam)
FILE=${IN_BAM[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}

# ----------------Commands------------------- #
# get number of mapped reads in millions
RPM=`grep -P ${BASE}'\t' library_sizes.txt | awk -F "\t" '{print $2}'`
RPM=`echo "scale=6; 1000000.0/"$RPM | bc`
BASE=${BASE}".scaled"

# bam to scaled bedGraph, bedGraph to bigWig
genomeCoverageBed -ibam $FILE -bg -scale $RPM -split -g $CHR_LENGTH > ${BASE}.bedGraph
wigToBigWig ${BASE}.bedGraph $CHR_LENGTH ${BASE}.bw
[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph
