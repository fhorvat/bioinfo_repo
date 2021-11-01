#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.06a.summarize_stats_and_tracks
#PBS -l select=ncpus=1:mem=1g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# read variables from source
source ./000.load_variables.sh

# script path
SCRIPT=04b.summarize_stats_and_tracks.script.R

# input for manual script
echo -e '\n'\
experiment=\'$EXPERIMENT\''\n'\
experiment_name=\'$EXPERIMENT_NAME\''\n'\
log_path=\'$LOG\''\n'

# ----------------Commands------------------- #
# run script
Rscript ${SCRIPT} \
--experiment ${EXPERIMENT} \
--experiment_name ${EXPERIMENT_NAME} \
--log_path ${LOG}
