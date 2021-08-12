#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=10g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.bam_to_psl
#PBS -J 0-4
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
INPUT_DIR=/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed/annotation
IN_SEQ=($(find $INPUT_DIR -name "*.bam"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/*/}
BASE=${BASE%.bam}

INPUT_DIR=../../mesAur
CHROM_SIZES=($(find ${INPUT_DIR} -name "*chrom.sizes"))

# ----------------Commands------------------- #
# bam to bigPsl
bamToPsl ${FILE} stdout | pslToBigPsl stdin stdout | LC_COLLATE=C sort -k1,1 -k2,2n > ${BASE}.bigPsl

# convert
bedToBigBed -extraIndex=name -as=/common/WORK/fhorvat/programi/UCSC/bigPsl.as -type=bed12+13 -tab ${BASE}.bigPsl ${CHROM_SIZES} ${BASE}.bb

# remove input
[ -f "${BASE}.bb" ] && rm ${BASE}.bigPsl

# set permission
chmod 744 ${BASE}.bb
