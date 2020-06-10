#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=2:mem=100g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.repeatMasker.genome
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script parameters
MEMORY=100G

# set input variables
INPUT_DIR=.
GENOME_FILE=(`ls ${INPUT_DIR}/*fa.gz`)
BASE=${GENOME_FILE#${INPUT_DIR}}
BASE=${BASE%.fa.gz}
OUT_DIR=.

# repeat masker parameters
RMSK_PAR="-e hmmer \
-s \
-species mouse \
-no_is \
-a \
-small \
-dir $OUT_DIR"

# ----------------Commands------------------- #
# unzip genome 
pigz -p 2 ${GENOME_FILE}

# run RepeatMasker
RepeatMasker $RMSK_PAR ${GENOME_FILE%.gz}

# zip genome
pigz -p 2 ${GENOME_FILE%.gz}
