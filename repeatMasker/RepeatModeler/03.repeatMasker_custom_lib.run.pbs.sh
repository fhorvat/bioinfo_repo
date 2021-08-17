#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=100g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03.repeatMasker.genome
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script parameters
MEMORY=100G
THREADS=12

# set input variables
INPUT_DIR=..
GENOME_FILE=($(find ${INPUT_DIR} -maxdepth 1 -name "*fasta"))
BASE=${GENOME_FILE#${INPUT_DIR}/}
BASE=${BASE%.fasta}

INPUT_DIR=.
FILE=($(find ${INPUT_DIR} -maxdepth 1 -name "*-families.fa"))
BASE=${FILE#${INPUT_DIR}}
BASE=${BASE%.fa}
OUT_DIR=./repeatMasker

mkdir $OUT_DIR

# repeat masker parameters
RMSK_PAR="-e rmblast \
-pa $THREADS \
-s \
-a \
-small \
-dir $OUT_DIR"

# ----------------Commands------------------- #
# run RepeatMasker
RepeatMasker $RMSK_PAR -lib ${FILE} ${GENOME_FILE}
