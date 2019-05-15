#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=10:mem=50g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.salmon_index
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script parameters
THREADS=10

# set input variables
IN_TRANS=`ls ../*.fasta`
OUT_PREFIX=salmon_index

# ----------------Commands------------------- #
# generate index
salmon index -t $IN_TRANS -i $OUT_PREFIX -k 31
