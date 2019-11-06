#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01a.find_element_coverage
#PBS -l select=ncpus=1:mem=20g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# number of threads
THREADS=1

# script path
SCRIPT=01b.find_element_coverage.R

# read variables from source
source ./00.load_variables.sh

# dataset path
DATASET_PATH=${DATASET_PATH}/4_merged_replicates

# features path
FEATURE_COORDINATES=${DOCUMENTATION_PATH}/LINE1_annotation.${L1_LIST}.GRanges.RDS

# input for manual script
echo -e experiment=\'$EXPERIMENT\''\n'\
single_end=$SINGLE_END'\n'\
threads=$THREADS'\n'\
dataset_path=\'$DATASET_PATH\''\n'\
documentation_path=\'$DOCUMENTATION_PATH \''\n'\
feature_coordinates=\'$FEATURE_COORDINATES\''\n'

# ----------------Commands------------------- #
# run script
Rscript $SCRIPT \
--experiment $EXPERIMENT \
--single_end $SINGLE_END \
--threads $THREADS \
--dataset_path $DATASET_PATH \
--documentation_path $DOCUMENTATION_PATH \
--feature_coordinates $FEATURE_COORDINATES
