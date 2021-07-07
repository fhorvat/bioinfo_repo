#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.04.bed_to_bigBed
#PBS -l select=ncpus=1:mem=1g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1

INPUT_DIR=.
IN_SEQ=(`find ${INPUT_DIR} -name "*.bed"`)
FILE=${IN_SEQ[0]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bed}

CHROM_SIZES=${INPUT_DIR}/*chrom.sizes

# ----------------Commands------------------- #
# bed to bigBed
bedToBigBed $FILE $CHROM_SIZES ${BASE}.bb
