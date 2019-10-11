#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01a.scale_merged_bigWig
#PBS -l select=ncpus=1:mem=10g
#PBS -J 0-2
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# set variables
SCRIPT=01b.scale_merged_bigWig.R

LOG=../../3_logs/log.read_stats.txt

INPUT_DIR=..
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 -name "*bw" -not -name "*scaled*"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bw}

# ----------------Commands------------------- #
# run script
Rscript $SCRIPT \
--bw_path $FILE \
--read_stats_path $LOG
