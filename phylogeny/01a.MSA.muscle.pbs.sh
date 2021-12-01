#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01a.MSA.muscle
#PBS -l select=ncpus=6:mem=10g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# input 
INPUT_DIR=..
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 \( -name "Dicer1.Metazoa.OrthoDB.20211117.taxonomy_and_sequences_PS_CDD.fasta" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fasta}

ALIGNER_NAME="muscle"
ALIGNER_SCRIPT=/common/WORK/kristian/bin/Muscle/muscle3.8.31_i86linux64

# ----------------Commands------------------- #
# MSA
${ALIGNER_SCRIPT} -in ${FILE} -out ${BASE}.${ALIGNER_NAME}.msa.fasta
