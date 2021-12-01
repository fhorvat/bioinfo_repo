#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01b.MSA.mafft
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

ALIGNER_NAME="mafft"
ALIGNER_SCRIPT=/common/WORK/kristian/bin/Miniconda/miniconda2/envs/qiime2-2019.4/bin/mafft

# ----------------Commands------------------- #
# MSA
${ALIGNER_SCRIPT} --maxiterate 1000 --genafpair --thread 6 ${FILE} > ${BASE}.${ALIGNER_NAME}.msa.fasta
