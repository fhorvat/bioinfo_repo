#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.gaps_in_assembly
#PBS -l select=ncpus=1:mem=20g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# set variables
SCRIPT=./01b.gaps_in_assembly.script.R

INPUT_DIR=../..
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 -name "*fa"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fa}

# input for manual script
echo -e '\n'\
fasta_path=\'$FILE\''\n'\
fasta_name=\'$BASE\''\n'

# ----------------Commands------------------- #
# run script
Rscript ${SCRIPT} --fasta_path ${FILE} --fasta_name ${BASE}
