#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.join_demultiplexed
#PBS -l select=ncpus=1:mem=20g
#PBS -j oe
#PBS -J 0-7
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# set input dir and find sample names
INPUT_DIR=.
IN_DIRS=($(find $INPUT_DIR -maxdepth 1 -type d -name "Undetermined*"))
DIR_NAMES=(${IN_DIRS[@]#${INPUT_DIR}/})
DIR_NAMES=(${DIR_NAMES[@]#Undetermined_from_})
DIR_NAMES=(${DIR_NAMES[@]%%_*})
DIR_NAMES=($(printf "%s\n" "${DIR_NAMES[@]}" | sort -u))

DATE=${DIR_NAMES[$PBS_ARRAY_INDEX]}

IN_SEQ=($(find $INPUT_DIR \( -name "*fastq.gz" \) | grep "$DATE" | sort | xargs echo ))

# just in case
printf "%s\n" "${IN_SEQ[@]}"

# ----------------Commands------------------- #
# join each demultiplexed samples to one file per sample per pairing
cat ${IN_SEQ[@]} > undetermined.${DATE}.fastq.gz
