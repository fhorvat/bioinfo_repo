#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.rmsk_expression
#PBS -l select=ncpus=1:mem=20g
#PBS -J 0-9
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# run parameters
THREADS=1

# script path
SCRIPT=01b.rmsk_expression.script.R

# input for script
INPUT_DIR=./bam_subset
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 -name "*.bam"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}

# input for manual script
echo -e '\n'\
bam_path=\'$FILE\''\n'\
bam_name=\'$BASE\''\n'

# ----------------Commands------------------- #
# run script
Rscript $SCRIPT \
--bam_path $FILE \
--bam_name $BASE
