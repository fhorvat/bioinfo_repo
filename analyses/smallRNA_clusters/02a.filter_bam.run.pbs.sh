#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02a.filter_bam
#PBS -l select=ncpus=12:mem=60g
#PBS -J 0-11
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script path
SCRIPT=02b.filter_bam.template.R

# read variables from source
source ./00.load_variables.sh

# override the number of threads
THREADS=12

# set path, get files, bam paths and names
IN_PATH=$MAPPED_PATH
IN_SEQ=($(find $IN_PATH -maxdepth 1 -name "*.bam"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE##${IN_PATH}/}
BASE=${BASE%.bam}

# get stats
STATS_PATH=./${BASE}.dt.RData

# set number of mismatches limit
NM_LIMIT=0

# input for manual script
echo -e '\n'\
experiment=\'$EXPERIMENT\''\n'\
threads=$THREADS'\n'\
genome_path=\'$GENOME_PATH\''\n'\
info_path=\'$INFO_PATH\''\n'\
exons_path=\'$EXONS_PATH\''\n'\
rmsk_path=\'$RMSK_PATH\''\n'\
mirbase_path=\'$MIRBASE_PATH\''\n'\
bam_path=\'$FILE\''\n'\
bam_name=\'$BASE\''\n'\
stats_path=\'$STATS_PATH\''\n'\
nm_limit=\'$NM_LIMIT\''\n'

# ----------------Commands------------------- #
# run script
Rscript $SCRIPT \
--experiment $EXPERIMENT \
--threads $THREADS \
--genome_path $GENOME_PATH \
--info_path $INFO_PATH \
--exons_path $EXONS_PATH \
--rmsk_path $RMSK_PATH \
--mirbase_path $MIRBASE_PATH \
--bam_path $FILE \
--bam_name $BASE \
--stats_path $STATS_PATH \
--nm_limit $NM_LIMIT
