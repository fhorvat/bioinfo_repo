#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.scale_bigWig
#PBS -l select=ncpus=1:mem=40g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
STAR_INDEX=/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1/STAR_index.2.7/sjdbOverhang_100
CHR_LENGTH=${STAR_INDEX}/chrNameLength.txt

INPUT_DIR=.
IN_SEQ=($(find ${INPUT_DIR} \( -name "*.bam" -and -name "*WT*" -not -name "*all_samples*" \)))
FILE=${IN_SEQ[0]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}

LOG=../3_logs/log.read_stats.txt

# ----------------Commands------------------- #
# just in case
awk -v PAT="$BASE" -F '\t' '$1 ~ PAT && $1 !~ "s_testis_Mov10l1_WT_8.5dpp_So820-M12_r3.SE" {print $1}' $LOG

# get number of mapped reads in millions
RPM=`awk -v PAT="$BASE" -F '\t' '$1 ~ PAT && $1 !~ "s_testis_Mov10l1_WT_8.5dpp_So820-M12_r3.SE" {sum += $8} END {print sum}' $LOG`
RPM=`echo "scale=6; 1000000.0/"$RPM | bc`
BASE=${BASE}.scaled

# bam to scaled bedGraph, bedGraph to bigWig
genomeCoverageBed -ibam $FILE -bg -scale $RPM -split -g $CHR_LENGTH > ${BASE}.bedGraph
wigToBigWig ${BASE}.bedGraph $CHR_LENGTH ${BASE}.bw
[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph
