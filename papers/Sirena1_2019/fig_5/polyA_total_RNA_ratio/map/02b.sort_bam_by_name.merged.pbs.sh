#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02b.sort_bam_by_name.merged
#PBS -l select=ncpus=10:mem=10g
#PBS -J 0-5
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
MEMORY=10GB
THREADS=10

INPUT_DIR=.
IN_SEQ=(`find $INPUT_DIR -name "*.genome.merged.bam"`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}

# ----------------Commands------------------- #
# sort by name
sambamba sort -m $MEMORY -o ${BASE}.sortedByName.bam -n -t $THREADS $FILE
