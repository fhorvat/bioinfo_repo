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

INPUT_DIR=/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 -name "*fasta"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fasta}

# ----------------Commands------------------- #
# run script
Rscript ${SCRIPT} --fasta_path ${FILE} 
