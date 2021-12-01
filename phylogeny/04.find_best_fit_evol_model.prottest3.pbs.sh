#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.04.find_best_fit_evol_model.prottest3
#PBS -l select=ncpus=6:mem=10g
#PBS -J 0-7
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# input 
INPUT_DIR=.
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 \( -name "*.msa.fasta-gb" -or -name "*.msa.fasta" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE}

SCRIPT="java -jar /common/WORK/fhorvat/programi/prottest3/prottest-3.4.2/prottest-3.4.2.jar"

# ----------------Commands------------------- #
# consensus tree
${SCRIPT} -i ${FILE} -o ${BASE}.AIC.best_fit.txt -all-distributions -AIC -threads 6
