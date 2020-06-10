#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=5g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N rmOutToGFF3
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
IN_DIR=.
OUT_DIR=.
FILE=($IN_DIR/*.fa.out)
BASE=${FILE#$IN_DIR/}
BASE=${BASE%.fa.out}
SCRIPT=/common/WORK/vfranke/bin/Repeatmasker/RepeatMasker/util/rmOutToGFF3.pl

# ----------------Commands------------------- #
$SCRIPT $FILE > ${OUT_DIR}/${BASE}".rmsk.gff"
