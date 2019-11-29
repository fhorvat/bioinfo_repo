#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.05a.0.class_reads
#PBS -l select=ncpus=1:mem=30g
#PBS -J 0-29
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# set path, get template, bams and names
IN_PATH=${PWD}
TEMPLATE=${IN_PATH}/05b.class_reads.template.R
SCRIPT_NAME=Rscript.${PBS_ARRAY_INDEX}.${TEMPLATE#${IN_PATH}/}

FILE=($IN_PATH/*.bam)
BAM=${FILE[$PBS_ARRAY_INDEX]}
BAM_NAME=`basename $BAM`
BAM_NAME=${BAM_NAME%.bam}
SCRIPT=${SCRIPT_NAME/template.R/${BAM_NAME}.R}

# ----------------Commands------------------- #
# change template, output to script, submit script
perl -pe "s|%BAM_PATH|$BAM|" $TEMPLATE > $SCRIPT
Rscript $SCRIPT

