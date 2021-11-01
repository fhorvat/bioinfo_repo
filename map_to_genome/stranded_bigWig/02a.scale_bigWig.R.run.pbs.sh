#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.scale_bigWig_in_R
#PBS -l select=ncpus=1:mem=40g
#PBS -J 0-11
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# set variables
SCRIPT=./02b.scale_bigWig.R.script.R

INPUT_DIR=.
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 -name "*bw" -not -name "*scaled*"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.plus.bw}
BASE=${BASE%.minus.bw}

INPUT_DIR_LOG=../3_logs
LOG=($(find $INPUT_DIR_LOG -name "log.read_stats.txt"))
SCALE_FACTOR=($(awk -F "\t" -vx="$BASE" '$1 == x {print $8 / 1000000}' $LOG))

# ----------------Commands------------------- #
# run script
Rscript $SCRIPT \
--bw_path $FILE \
--scale_factor $SCALE_FACTOR
