#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=1g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.split_fasta_by_reference
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
IN_DIR=.
IN_SEQ=($(find ${IN_DIR} -name "*fasta"))
FILE=${IN_SEQ[0]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.fasta}.split

# create dir
mkdir $BASE

# ----------------Commands------------------- #
# split fasta on individual references
faSplit byname $FILE ${BASE}/
