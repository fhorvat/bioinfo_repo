#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01b.sort_bam_by_name.genome
#PBS -l select=ncpus=10:mem=10g
#PBS -J 0-%N_SAMPLES
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
MEMORY=10GB
THREADS=10

INPUT_DIR=.
IN_SEQ=(`find $INPUT_DIR -name "*.genome.bam"`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}

# ----------------Commands------------------- #
# sort by name
sambamba sort -m $MEMORY -o ${BASE}.sortedByName.bam -n -t $THREADS $FILE
