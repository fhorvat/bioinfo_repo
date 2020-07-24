#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.filter_unique
#PBS -l select=ncpus=6:mem=30g
#PBS -j oe
#PBS -J 0-8
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6

# files
INPUT_DIR=.
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 -name "*.bam" -not -name "*unique.bam"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}.unique

# ----------------Commands------------------- #
# remove unmapped reads, duplicates and reads with multiple alignments
sambamba view -h -t $THREADS -f bam \
-F "[XS] == null and not unmapped and not duplicate" \
$FILE > ${BASE}.bam
