#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.annotate_assembly
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6
MEMORY=20G

INPUT_DIR=.
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 -name "*fasta" -not -name "*l1000*"))
FILE=${IN_SEQ[0]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fasta}

SCRIPT=/common/WORK/fhorvat/Projekti/Svoboda/scripts/other_scripts/filterFastaByLength.pl

# ----------------Commands------------------- #
# filter fasta by length
perl $SCRIPT 1000 ${FILE} > ${BASE}.l1000.fasta

# find most likely ORF candidates
TransDecoder.LongOrfs -t ${BASE}.l1000.fasta 
