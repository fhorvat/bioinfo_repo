#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=50g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.RepeatModeler
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script parameters
MEMORY=50G
THREADS=6

# set input variables
INPUT_DIR=..
GENOME_FILE=($(find ${INPUT_DIR} -maxdepth 1 -name "*fasta"))
BASE=${GENOME_FILE#${INPUT_DIR}/}
BASE=${BASE%.fasta}

# set script path
SCRIPT=/common/WORK/fhorvat/programi/RepeatModeler/RepeatModeler-2.0.1/RepeatModeler

# ----------------Commands------------------- #
# create RepeatModeler database
$SCRIPT -database ${BASE} -pa $THREADS -LTRStruct >& ${BASE}.out
