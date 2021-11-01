#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=50g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.bsbolt_index
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script parameters
THREADS=6

# set input variables
INPUT_DIR=.
IN_SEQ=$(find $INPUT_DIR -maxdepth 1 -name "hamster.sequel.draft-20200302.arrow.fasta")
FILE=${IN_SEQ[0]}
BASE=${BASE#${INPUT_DIR}/}
BASE=${BASE%.fasta}
BASE=${BASE%.fasta}

# ----------------Commands------------------- #
# generate index
bsbolt Index -G $FILE -DB .
