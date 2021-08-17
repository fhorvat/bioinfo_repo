#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=100g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.repeatMasker.genome
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script parameters
MEMORY=100G
THREADS=12

# set input variables
INPUT_DIR=.
GENOME_FILE=($(find ${INPUT_DIR} -name "*fasta"))
BASE=${GENOME_FILE#${INPUT_DIR}}
BASE=${BASE%.fasta}
OUT_DIR=./repeatMasker

mkdir $OUT_DIR

# repeat masker parameters
RMSK_PAR="-e hmmer \
-pa $THREADS \
-s \
-species mouse \
-no_is \
-a \
-small \
-dir $OUT_DIR"

# ----------------Commands------------------- #
# run RepeatMasker
RepeatMasker $RMSK_PAR $GENOME_FILE
