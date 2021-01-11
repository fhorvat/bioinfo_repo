#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.get_list_of_converted_reads
#PBS -l select=ncpus=1:mem=10g
#PBS -j oe
#PBS -J 0-1
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script path
SCRIPT=01b.get_list_of_converted_reads.chunked.script.R

# get fastq
IN_DIR=.
IN_SEQ=($(find ${IN_DIR} \( -name "s_mouse*fastq" -not -name "*converted.fastq" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.fastq}

# input for manual script
echo -e '\n'\
fastq_path=\'$FILE\''\n'

# ----------------Commands------------------- #
# run script
Rscript $SCRIPT \
--fastq_path $FILE
