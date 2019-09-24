#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=10:mem=50g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.repeatMasker.rodents
#PBS -J 2-3
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script parameters
THREADS=10
MEMORY=50G

# set input variables
INPUT_DIR=..
GENOME_FILES=(${INPUT_DIR}/*.fa)
GENOME_FILE=${GENOME_FILES[$PBS_ARRAY_INDEX]}
BASE=${GENOME_FILE#${INPUT_DIR}}
BASE=${BASE%.fa}
OUT_DIR=./${BASE}

# create out dir
mkdir -p $OUT_DIR

# repeat masker parameters
RMSK_PAR="-pa $THREADS \
-e hmmer \
-s \
-species mouse \
-no_is \
-a \
-small \
-dir $OUT_DIR"

# ----------------Commands------------------- #
# run RepeatMasker
RepeatMasker $RMSK_PAR $GENOME_FILE
