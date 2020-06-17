#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.00.bwa_index
#PBS -l select=ncpus=1:mem=10g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
INPUT_DIR=$PWD
IN_FASTA=${INPUT_DIR}/L1.Mus.20200507.dfam.fa

# ----------------Commands------------------- #
# BWA
bwa index $IN_FASTA
