#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=5g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.05.assembly_stats
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1
MEMORY=5G

IN_DIR_TRANS=.
IN_SEQ_TRANS=(${IN_DIR_TRANS}/*.fasta)
FILE_TRANS=${IN_SEQ_TRANS[0]}

BASE=${FILE_TRANS#${IN_DIR_TRANS}/}
BASE=${BASE}.stats

SCRIPT=/common/WORK/fhorvat/programi/Trinity/trinityrnaseq-Trinity-v2.8.4/util/TrinityStats.pl

# ----------------Commands------------------- #
$SCRIPT $FILE_TRANS > $BASE
