#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.04.assemble_contigs.wtdbg2
#PBS -l select=ncpus=1:mem=5g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1

INPUT_DIR=.
IN_SEQ=(`find $INPUT_DIR -maxdepth 1 -name "*fasta" -not -name "*assembled*"`)
FILE=${IN_SEQ[0]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fasta}

OUT_DIR=./${BASE}
SCRIPT=/common/WORK/fhorvat/programi/wtdbg2/wtdbg2/wtdbg2
SCRIPT_2=/common/WORK/fhorvat/programi/wtdbg2/wtdbg2/wtpoa-cns
# ----------------Commands------------------- #
# assemble
$SCRIPT \
-k 23 -p 0 -K 1000 -A -S 1 -s 0.05 -X 1 -e 1 -L 1 \
-l 100 -m 100 \
-g 10000k \
-t $THREADS \
-fo assembly.${BASE} \
-i $FILE

# derive consensus
$SCRIPT_2 \
-t $THREADS \
-i assembly.${BASE}.ctg.lay.gz \
-fo assembly.${BASE}.raw.fa
