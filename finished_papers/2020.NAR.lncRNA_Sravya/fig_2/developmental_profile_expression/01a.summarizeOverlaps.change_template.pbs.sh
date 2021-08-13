#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01a.summarizeOverlaps.PE
#PBS -l select=ncpus=10:mem=100g
#PBS -J 0-1
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# set path, get template, bams and names
IN_PATH=.
TEMPLATE=`ls ${IN_PATH}/*.summarizeOverlaps.template.R`
SCRIPT_NAME=Rscript.${TEMPLATE#${IN_PATH}/}

FILES=(Fugaku \
Veselovska_2015_GenomeBiol_GSE70116)

EXPERIMENT=${FILES[$PBS_ARRAY_INDEX]}
SCRIPT=${SCRIPT_NAME/template.R/${EXPERIMENT}.R}

# ----------------Commands------------------- #
# change template, output to script
perl -pe "s|%EXPERIMENT|$EXPERIMENT|" $TEMPLATE > $SCRIPT

# run script
Rscript $SCRIPT
