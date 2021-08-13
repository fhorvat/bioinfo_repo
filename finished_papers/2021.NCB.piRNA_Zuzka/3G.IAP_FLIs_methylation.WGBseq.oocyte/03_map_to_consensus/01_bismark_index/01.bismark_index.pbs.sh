#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=24:mem=50g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.bismark_index
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script parameters
THREADS=12

# set input variables
INPUT_DIR=.
IN_SEQ=$(find $INPUT_DIR -name "*.fasta" -and -name "IAP.potentially_young.ordered_by_ORFs.20201031.IAPLTR3.consensus.fasta")
FILE=${IN_SEQ[0]}
BASE=${BASE#${INPUT_DIR}/}
BASE=${BASE%.fasta}

# ----------------Commands------------------- #
# generate index
bismark_genome_preparation --verbose --bowtie2 --parallel ${THREADS} ${INPUT_DIR}
