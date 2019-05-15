#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.04a.class_reads
#PBS -l select=ncpus=1:mem=20g
#PBS -J 0-%N_SAMPLES
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# set path, get template, bams and names
TEMPLATE=(04b.class_reads.template.R)
SCRIPT_NAME=Rscript.${TEMPLATE}

IN_PATH=.
FILE=(`ls ${IN_PATH}/*.sortedByName.bam`)
BAM=${FILE[$PBS_ARRAY_INDEX]}
BAM_NAME=${BAM#${IN_PATH}/}
BAM_NAME=${BAM_NAME%.sortedByName.bam}
SCRIPT=${SCRIPT_NAME/template.R/${BAM_NAME}.R}

# ----------------Commands------------------- #
# change template, output to script, submit script
perl -pe "s|%BAM_PATH|$BAM|" $TEMPLATE > $SCRIPT
Rscript $SCRIPT
