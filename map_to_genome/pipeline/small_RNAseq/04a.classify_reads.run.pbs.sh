#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.04a.classify_reads
#PBS -l select=ncpus=1:mem=5g
#PBS -j oe
#PBS -J 0-%N_TOTAL
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script path
SCRIPT=04b.classify_reads.script.R

# read variables from source
source ./000.load_variables.sh

# files path
IN_PATH=.
FILE=($(find ${IN_PATH} -maxdepth 1 -name "*.bam"))
BAM_PATH=${FILE[$PBS_ARRAY_INDEX]}
BAM_NAME=${BAM_PATH#${IN_PATH}/}
BAM_NAME=${BAM_NAME%.bam}

# input for manual script
echo -e '\n'\
single_end=\'$SINGLE_END\''\n'\
threads=\'$THREADS\''\n'\
bam_path=\'$BAM_PATH\''\n'\
bam_name=\'$BAM_NAME\''\n'\
features_exons=\'$FEATURES_EXONS\''\n'\
features_rmsk=\'$FEATURES_RMSK\''\n'\
features_geneInfo=\'$FEATURES_GENEINFO\''\n'\
features_mirbase=\'$FEATURES_MIRBASE\''\n'

# ----------------Commands------------------- #
# run script
Rscript $SCRIPT \
--single_end $SINGLE_END \
--threads $THREADS \
--bam_path $BAM_PATH \
--bam_name $BAM_NAME \
--features_exons $FEATURES_EXONS \
--features_geneInfo $FEATURES_GENEINFO \
--features_rmsk $FEATURES_RMSK \
--features_mirbase $FEATURES_MIRBASE
