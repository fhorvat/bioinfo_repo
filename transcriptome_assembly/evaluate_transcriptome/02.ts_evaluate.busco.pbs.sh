#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=20:mem=30g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.busco_ts_evaluate
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script options
THREADS=20
BUSCO=/common/WORK/fhorvat/programi/busco/busco-master/scripts/run_BUSCO.py
BUSCO_LINEAGE=/common/WORK/fhorvat/programi/busco/busco-master/lineages/mammalia_odb9
BUSCO_BASE=`basename $BUSCO_LINEAGE`
BUSCO_BASE=${BUSCO_BASE%_odb9}

# fasta files
INPUT_DIR=.
IN_SEQ=(${INPUT_DIR}/*.fasta)
FILE=${IN_SEQ[0]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fasta}.busco

# ----------------Commands------------------- #
# evaluate fasta files
python3 $BUSCO \
-i $FILE \
-o $BASE.$BUSCO_BASE \
-l $BUSCO_LINEAGE \
-m tran \
-c $THREADS
