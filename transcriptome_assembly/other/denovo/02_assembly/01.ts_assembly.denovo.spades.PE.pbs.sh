#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=250g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.spades_ts_assembly
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=250G

INPUT_DIR=..
IN_SEQ=(`ls ${INPUT_DIR}/*.txt.gz`)
FILES=$(printf ",%s" ${IN_SEQ[@]})
FILES=${FILES:1}

BASE="pwd"
OUT_DIR=${BASE}_trinity_assembly

SPADES=/common/WORK/dglavas/programs/SPAdes-3.13.0-Linux/bin/rnaspades.py

SPADES_PAR="-o . \
-t $THREADS \
-m $MEMORY \
-k 21,23"
 
# ----------------Commands------------------- #
# run transcriptome assembly
$SPADES \
--pe1-1 $IN_FOLDER/Esu1RNA/Esu1RNA_1_adtrimmed.fq.gz \
--pe1-2 $IN_FOLDER/Esu1RNA/Esu1RNA_2_adtrimmed.fq.gz \
--pe2-1 $IN_FOLDER/Esu2RNA/Esu2RNA_1_adtrimmed.fq.gz \
--pe2-2 $IN_FOLDER/Esu2RNA/Esu2RNA_2_adtrimmed.fq.gz \
$SPADES_PAR | awk '{ print strftime("%Y-%m-%d %H:%M:%S"), $0; fflush(); }' &> spades_timestamp.log

