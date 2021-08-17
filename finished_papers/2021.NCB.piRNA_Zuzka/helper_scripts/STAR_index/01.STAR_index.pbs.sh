#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=18:mem=200g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.STAR_index
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script parameters
THREADS=18

# set input variables
IN_DIR=..
IN_GENOME=($(find $IN_DIR -name "*.fasta"))
OUTDIR=.

# create outdir
mkdir $OUTDIR 2> /dev/null

# ----------------Commands------------------- #
# generate index
STAR.2.7 \
--runThreadN $THREADS \
--runMode genomeGenerate \
--genomeDir $OUTDIR \
--genomeFastaFiles $IN_GENOME \
--limitGenomeGenerateRAM=50000000000
