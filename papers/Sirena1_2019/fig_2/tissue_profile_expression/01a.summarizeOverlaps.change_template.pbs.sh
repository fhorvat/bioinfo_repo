#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01a.summarizeOverlaps.PE
#PBS -l select=ncpus=10:mem=100g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# set path, get template, bams and names
IN_PATH=.
TEMPLATE=`ls ${IN_PATH}/*.summarizeOverlaps.template.R`
SCRIPT_NAME=Rscript.${TEMPLATE#${IN_PATH}/}

FILES=(ENCODE_2014_Nature_GSE49417)

EXPERIMENT=${FILES[0]}
SCRIPT=${SCRIPT_NAME/template.R/${EXPERIMENT}.R}

# ----------------Commands------------------- #
# change template, output to script
perl -pe "s|%EXPERIMENT|$EXPERIMENT|" $TEMPLATE > $SCRIPT

# run script
Rscript $SCRIPT
