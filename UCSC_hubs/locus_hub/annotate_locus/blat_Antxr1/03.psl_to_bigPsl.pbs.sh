#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03.psl_to_bigPsl
#PBS -l select=ncpus=1:mem=1g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1

INPUT_DIR=.
IN_SEQ=($(find ${INPUT_DIR} \( -name "*.psl" -not -name "*score.psl" -not -name "*all*" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.psl}

CHROM_SIZES_INPUT_DIR=../chrom_sizes
CHROM_SIZES=(${CHROM_SIZES_INPUT_DIR}/${BASE/.*/*chrom.sizes})

# ----------------Commands------------------- #
# prepare psl for bigPsl conversion (sort)
pslToBigPsl $FILE stdout | sort -k1,1 -k2,2n > ${BASE}.bigPslInput

# convert
bedToBigBed -as=/common/WORK/fhorvat/programi/UCSC/bigPsl.as -type=bed12+13 -tab ${BASE}.bigPslInput ${CHROM_SIZES} ${BASE}.bb

# remove input
[ -f "${BASE}.bb" ] && rm ${BASE}.bigPslInput

