#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03a.class_reads.change_template
#PBS -l select=ncpus=1:mem=100g
#PBS -J 0-%N_SAMPLES
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# set path, get template, bams and names
IN_PATH=.
TEMPLATE=`ls ${IN_PATH}/*.class_reads.template.R`
SCRIPT_NAME=Rscript.${TEMPLATE#${IN_PATH}/}

FILE=(`ls ${IN_PATH}/*genome.Aligned.sortedByCoord.out.bam`)
BAM=${FILE[$PBS_ARRAY_INDEX]}
BAM_NAME=${BAM#${IN_PATH}/}
BAM_NAME=${BAM_NAME%.genome.Aligned.sortedByCoord.out.bam}
SCRIPT=${SCRIPT_NAME/template.R/${BAM_NAME}.R}

# ----------------Commands------------------- #
# change template, output to script, submit script
perl -pe "s|%BAM_PATH|$BAM|" $TEMPLATE > $SCRIPT
Rscript $SCRIPT
