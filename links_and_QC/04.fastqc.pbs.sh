#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=1g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03.fastqc
#PBS -j oe
#PBS -J 0-15
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
INPUT_DIR=../Links
OUTPUT_DIR=.
IN_SEQ=(`find $INPUT_DIR -regex ".*[txt,fastq].gz"`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}

# ----------------Commands------------------- #
fastqc --outdir $OUTPUT_DIR $FILE
