#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01a.summarizeOverlaps.PE
#PBS -l select=ncpus=10:mem=100g
##PBS -J 0-5
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# set path, get template, bams and names
IN_PATH=.
TEMPLATE=`ls ${IN_PATH}/*.summarizeOverlaps.template.R`
SCRIPT_NAME=Rscript.${TEMPLATE#${IN_PATH}/}

FILES=(Stein_2015_PLoSGenet_GSE57514)

EXPERIMENT=${FILES[0]}
SCRIPT=${SCRIPT_NAME/template.R/${EXPERIMENT}.R}
OUT_PATH=`echo $PWD`

# ----------------Commands------------------- #
# change template, output to script
perl -pe "s|%EXPERIMENT|$EXPERIMENT|" $TEMPLATE > $SCRIPT
perl -pi -e "s|%OUT_PATH|$OUT_PATH|" $SCRIPT

# run script
Rscript $SCRIPT
