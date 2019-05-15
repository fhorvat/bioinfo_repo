#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.05.sum_reads
#PBS -l select=ncpus=1:mem=1g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script path
SCRIPT=05b.sum_reads.template.R

# ----------------Commands------------------- #
# run script
Rscript $SCRIPT
