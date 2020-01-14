#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.FPM_tables
#PBS -l select=ncpus=1:mem=20g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script path
SCRIPT=02b.FPM_tables.script.R

# read variables from source
source ./00.load_variables.sh

# counts path
COUNTS_PATH=($(find . -maxdepth 1 -name "*.counts.txt"))

# input for manual script
echo -e '\n'\
experiment=\'$EXPERIMENT\''\n'\
threads=$THREADS'\n'\
mapped_path=\'$MAPPED_PATH\''\n'\
documentation_path=\'$DOCUMENTATION_PATH\''\n'\
features_coordinates=\'$FEATURES_COORDINATES\''\n'\
features_name=\'$FEATURES_NAME\''\n'\
genes_info_path=\'$GENES_INFO_PATH\''\n'\
grouping_variables=\'$GROUPING_VARIABLES\''\n'\
counts_path=\'$COUNTS_PATH\''\n'

# ----------------Commands------------------- #
# run script
Rscript $SCRIPT \
--experiment $EXPERIMENT \
--threads $THREADS \
--mapped_path $MAPPED_PATH \
--documentation_path $DOCUMENTATION_PATH \
--features_coordinates $FEATURES_COORDINATES \
--features_name $FEATURES_NAME \
--genes_info_path $GENES_INFO_PATH \
--grouping_variables $GROUPING_VARIABLES \
--counts_path $COUNTS_PATH
