#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=10g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.soft_mask_genome
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script parameters
THREADS=1

# set input variables
IN_DIR=..
IN_GENOME=($(find $IN_DIR -maxdepth 1 -name "*.fasta"))
BASE=${IN_GENOME#${IN_DIR}/}
BASE=${BASE%.fasta}

IN_DIR=.
IN_BED=($(find $IN_DIR -maxdepth 1 -name "*.bed"))

# ----------------Commands------------------- #
# soft mask
bedtools maskfasta -soft -fi ${IN_GENOME} -bed ${IN_BED} -fo ${BASE}.sm.fasta
