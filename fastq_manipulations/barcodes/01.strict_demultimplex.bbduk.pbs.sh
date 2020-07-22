#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.strict_demultiplex
#PBS -l select=ncpus=6:mem=30g
#PBS -J 0-9
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6
MEMORY=30g

INPUT_DIR=..
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 -name "*.txt.gz"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.txt.gz}

# get barcodes
SAMPLE_TABLE=$(find ../../../Documentation -name "*sampleTable.csv")
BARCODE=$(awk -F "," -vx="$BASE" '$1 == x {print $7}' $SAMPLE_TABLE)

# ----------------Commands------------------- #
bbduk.sh in=${FILE} out=${BASE}.strict.txt.gz barcodes=${BARCODE} barcodefilter chastityfilter threads=$THREADS
