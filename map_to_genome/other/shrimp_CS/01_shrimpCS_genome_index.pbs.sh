#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N shrimpCS
#PBS -l select=ncpus=1:mem=50g
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
REF_DIR=/common/WORK/fhorvat/reference/pig/susScr3
SCRIPT=/common/WORK/fhorvat/programi/SHRiMP2/SHRiMP_2_2_3/utils/project-db.py

# ----------------Commands------------------- #
# create projections of genome
$SCRIPT --shrimp-mode cs ${REF_DIR}/UCSC/susScr3.fa
