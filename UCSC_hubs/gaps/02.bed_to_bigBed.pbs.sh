#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.bed_to_bigBed
#PBS -l select=ncpus=1:mem=1g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1

INPUT_DIR=.
IN_SEQ=($(find $INPUT_DIR -name "*.bed" ! -name "*sorted.bed"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bed}

INPUT_DIR=../../mesAur
CHROM_SIZES=($(find ${INPUT_DIR} -name "*chrom.sizes"))

# ----------------Commands------------------- #
# prepare psl for bigPsl conversion (sort)
LC_COLLATE=C sort -k1,1 -k2,2n $FILE > ${BASE}.sorted.bed

# convert
bedToBigBed ${BASE}.sorted.bed ${CHROM_SIZES} ${BASE}.bb

# remove input
[ -f "${BASE}.bb" ] && rm ${BASE}.sorted.bed
#[ -f "${BASE}.bb" ] && rm ${BASE}.bed

# set permission
chmod 744 ${BASE}.bb
