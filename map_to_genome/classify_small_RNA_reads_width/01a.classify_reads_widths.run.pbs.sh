#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01a.classify_reads_widths
#PBS -l select=ncpus=1:mem=50g
#PBS -j oe
#PBS -J 0-2
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
features_exons=\'$FEATURES_EXONS\''\n'\
features_rmsk=\'$FEATURES_RMSK\''\n'\
features_geneInfo=\'$FEATURES_GENEINFO\''\n'\
features_mirbase=\'$FEATURES_MIRBASE\''\n'

# ----------------Commands------------------- #
# run script
Rscript $SCRIPT \
--single_end $SINGLE_END \
--bam_path $BAM_PATH \
--bam_name $BAM_NAME \
--features_exons $FEATURES_EXONS \
--features_geneInfo $FEATURES_GENEINFO \
--features_rmsk $FEATURES_RMSK \
--features_mirbase $FEATURES_MIRBASE
