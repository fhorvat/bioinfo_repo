#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.bed_to_bigBed
#PBS -l select=ncpus=1:mem=1g
#PBS -J 0-5
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1

INPUT_DIR=.
IN_SEQ=(`find ${INPUT_DIR} -name "*.bed"`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bed}

CHROM_SIZES_INPUT_DIR=../chrom_sizes
CHROM_SIZES=(${CHROM_SIZES_INPUT_DIR}/${BASE/.*/*chrom.sizes})

# ----------------Commands------------------- #
# sort
sort -k1,1 -k2,2n $FILE > ${BASE}.bigBedInput

# convert
bedToBigBed ${BASE}.bigBedInput ${CHROM_SIZES} ${BASE}.bb

# remove input
[ -f "${BASE}.bb" ] && rm ${BASE}.bigBedInput

