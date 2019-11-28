#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=40g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.bbmap_index
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script parameters
MEMORY=40G
THREADS=12

# set input variables
IN_GENOME=../mm10.fa.gz
OUTDIR=.

# ----------------Commands------------------- #
# generate index
bbmap.sh -Xmx$MEMORY threads=$THREADS ref=$IN_GENOME

