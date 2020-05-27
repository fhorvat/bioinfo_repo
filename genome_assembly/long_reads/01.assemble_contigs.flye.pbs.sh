#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.assemble_contigs.flye
#PBS -l select=ncpus=1:mem=5g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1

INPUT_DIR=.
IN_SEQ=(`find $INPUT_DIR -maxdepth 1 -name "*fasta"`)
FILE=${IN_SEQ[0]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fasta}

OUT_DIR=./${BASE}
SCRIPT=/common/WORK/kristian/bin/Flye/Flye-2.5/bin/flye

# ----------------Commands------------------- #
# assemble
$SCRIPT --subassemblies $FILE --genome-size 1m --out-dir $OUT_DIR --threads $THREADS --min-overlap 500
