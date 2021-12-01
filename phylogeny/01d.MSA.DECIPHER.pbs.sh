#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01d.MSA.DECIPHER
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

ALIGNER_NAME="DECIPHER"
ALIGNER_SCRIPT="Rscript /common/WORK/fhorvat/Projekti/Svoboda/scripts/R_scripts/DECIPHER.R"

# ----------------Commands------------------- #
# MSA
${ALIGNER_SCRIPT} --in_file ${FILE} --out_file ${BASE}.${ALIGNER_NAME}.msa.fasta --threads 6 
