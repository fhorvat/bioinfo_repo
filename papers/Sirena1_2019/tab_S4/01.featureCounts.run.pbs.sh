#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01c.featureCounts
#PBS -l select=ncpus=12:mem=40g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# read variables from source
source ./00.load_variables.sh

# override the number of threads
THREADS=12

# files path
IN_SEQ=($(find $MAPPED_PATH -maxdepth 1 -name "*.bam"))
IN_SEQ=`printf "%s " ${IN_SEQ[@]}`

# script parameters
SCRIPT_PARAMS="-T $THREADS \
--fracOverlap 1 \
-p \
-F SAF \
-M \
-O"

# ----------------Commands------------------- #
# run script
featureCounts -a $FEATURES_COORDINATES -o ${FEATURES_NAME}.counts.txt $IN_SEQ $SCRIPT_PARAMS
