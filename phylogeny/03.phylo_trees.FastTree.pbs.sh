#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03.phylo_trees.FastTree
#PBS -l select=ncpus=6:mem=1g
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

SCRIPT=/common/WORK/fhorvat/programi/FastTree/FastTreeMP

# ----------------Commands------------------- #
# maximum likelihood phylogenetic trees
${SCRIPT} ${FILE} > ${BASE}.tree
