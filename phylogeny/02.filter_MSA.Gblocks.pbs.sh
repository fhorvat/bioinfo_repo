#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.filter_MSA.Gblocks
#PBS -l select=ncpus=1:mem=1g
#PBS -J 0-3
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# input 
INPUT_DIR=.
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 \( -name "*.msa.fasta" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fasta}

SCRIPT=/common/WORK/fhorvat/programi/Gblocks/Gblocks_0.91b/Gblocks

# ----------------Commands------------------- #
# get number of sequences
SEQ_COUNT=$(grep -c ">" ${FILE})
SEQ_COUNT=$(echo "scale=0; 0.6*"${SEQ_COUNT}" / 1" | bc)

# clean MSA
${SCRIPT} ${FILE} -t=p -b2=${SEQ_COUNT} -b3=10 -b4=5 -b5=a -p=n
