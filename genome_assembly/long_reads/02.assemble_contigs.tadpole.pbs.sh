#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.assemble_contigs.tadpole
#PBS -l select=ncpus=1:mem=40g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1
MEMORY=40G

INPUT_DIR=.
IN_SEQ=(`find $INPUT_DIR -maxdepth 1 \( -name "*fasta" -not -name "*assembled*" \) `)
FILE=${IN_SEQ[0]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fasta}

# ----------------Commands------------------- #
# assemble
tadpole.sh -Xmx$MEMORY threads=$THREADS in=$FILE out=${BASE}.assembled.fasta overwrite=true -k=150
