#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.psl_to_bigPsl
#PBS -l select=ncpus=1:mem=1g
#PBS -j oe
#PBS -J 0-1
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1

INPUT_DIR=.
IN_SEQ=($(find ${INPUT_DIR} \( -name "*.psl" -not -name "*score.psl" -not -name "*all*" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.psl}

CHROM_SIZES_INPUT_DIR=/common/DB/genome_reference/Mollusca/Arion_vulgaris.Schrodl/STAR_index.2.7/sjdbOverhang_100
CHROM_SIZES=(${CHROM_SIZES_INPUT_DIR}/chrNameLength.txt)

# ----------------Commands------------------- #
# prepare psl for bigPsl conversion (sort)
pslToBigPsl $FILE stdout | LC_COLLATE=C sort -k1,1 -k2,2n > ${BASE}.bigPslInput

# convert
bedToBigBed -as=/common/WORK/fhorvat/programi/UCSC/bigPsl.as -type=bed12+13 -tab ${BASE}.bigPslInput ${CHROM_SIZES} ${BASE}.bb

# remove input
[ -f "${BASE}.bb" ] && rm ${BASE}.bigPslInput

