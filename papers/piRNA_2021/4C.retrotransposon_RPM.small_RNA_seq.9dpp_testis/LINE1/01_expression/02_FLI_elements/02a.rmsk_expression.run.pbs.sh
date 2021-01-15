#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.rmsk_expression
#PBS -l select=ncpus=1:mem=10g
#PBS -J 0-4
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# run parameters
THREADS=1

# script path
SCRIPT=02b.rmsk_expression.script.R

# input for script
INPUT_DIR=./bam_subset
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 -name "*.bam"))
FILE_BAM=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE_BAM=${FILE_BAM#${INPUT_DIR}/}
BASE_BAM=${BASE_BAM%.bam}

# bed input
INPUT_DIR=../..
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 -name "*.FLI_elements.bed"))
FILE_BED=${IN_SEQ[0]}
BASE_BED=${FILE_BED#${INPUT_DIR}/}
BASE_BED=${BASE_BED%.bed}

# input for manual script
echo -e '\n'\
bam_path=\'$FILE_BAM\''\n'\
bam_name=\'$BASE_BAM\''\n'\
bed_path=\'$FILE_BED\''\n'

# ----------------Commands------------------- #
# run script
Rscript $SCRIPT \
--bam_path $FILE_BAM \
--bam_name $BASE_BAM \
--bed_path $FILE_BED
