#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01a.classify_reads_widths
#PBS -l select=ncpus=1:mem=50g
#PBS -j oe
#PBS -J 0-4
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script path
SCRIPT=01b.classify_reads_widths.script.R

# read variables from source
source ./000.load_variables.sh

# files path
BAM_PATH=${IN_SEQ[$PBS_ARRAY_INDEX]}
BAM_NAME=${BAM_PATH#${INPUT_DIR}/}
BAM_NAME=${BAM_NAME%.bam}

# input for manual script
echo -e '\n'\
single_end=\'$SINGLE_END\''\n'\
bam_path=\'$BAM_PATH\''\n'\
bam_name=\'$BAM_NAME\''\n'\
experiment_name=\'$EXPERIMENT_NAME\''\n'\
features_ensembl=\'$FEATURES_ENSEMBL\''\n'\
features_rmsk=\'$FEATURES_RMSK\''\n'

# ----------------Commands------------------- #
# run script
Rscript $SCRIPT \
--single_end $SINGLE_END \
--bam_path $BAM_PATH \
--bam_name $BAM_NAME \
--experiment_name $EXPERIMENT_NAME \
--features_ensembl $FEATURES_ENSEMBL \
--features_rmsk $FEATURES_RMSK
