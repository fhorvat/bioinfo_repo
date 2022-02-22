#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.featureCounts
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

# script parameters - single or paired end
if [ $SINGLE_END == TRUE ]
then
   SCRIPT_PARAMS="-T $THREADS -F GTF -M -O --fraction"
else
   SCRIPT_PARAMS="-T $THREADS -F GTF -M -O --fraction -p"
fi

# ----------------Commands------------------- #
# run script
echo "Running feature counts with following parameters:" $SCRIPT_PARAMS
echo "Single end reads are" $SINGLE_END
featureCounts -a $FEATURES_COORDINATES -o ${FEATURES_NAME}.counts.txt $IN_SEQ $SCRIPT_PARAMS
