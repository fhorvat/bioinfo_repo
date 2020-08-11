#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.featureCounts
#PBS -l select=ncpus=12:mem=40g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# number of threads
THREADS=12

# read variables from source
source ./00.load_variables.sh

# feature coordinates 
FEATURE_COORDINATES=${DOCUMENTATION_PATH}/LINE1_annotation.${L1_LIST}.saf

# files path
IN_SEQ=($(find $DATASET_PATH -maxdepth 1 -name "*.bam"))
IN_SEQ=`printf "%s " ${IN_SEQ[@]}`

# script parameters
SCRIPT_PARAMS="-T $THREADS \
-F SAF \
-M \
-p \
-O"

# ----------------Commands------------------- #
# run script
featureCounts -a $FEATURE_COORDINATES -o LINE1_annotation.${L1_LIST}.${EXPERIMENT}.counts.txt $IN_SEQ $SCRIPT_PARAMS
