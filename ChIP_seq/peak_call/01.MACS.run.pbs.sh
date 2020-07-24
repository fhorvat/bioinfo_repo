#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.MACS
#PBS -l select=ncpus=1:mem=40g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# read variables from source
source ./00.load_variables.sh

# override the number of threads
THREADS=1

# outdir
OUTDIR=.

# treated in
IN_SEQ=($(find $MAPPED_PATH -maxdepth 1 -name "*unique.bam" -and -name "*ChIP*" -not -name "*input*"))
IN_SEQ=`printf "%s " ${IN_SEQ[@]}`

# input in
IN_INPUT=($(find $MAPPED_PATH -maxdepth 1 -name "*unique.bam" -and -name "*ChIP*" -and -name "*input*"))
IN_INPUT=`printf "%s " ${IN_INPUT[@]}`

# ----------------Commands------------------- #
# run script
macs2 callpeak -t ${IN_SEQ} -c ${IN_INPUT} -n ${EXPERIMENT}.ChIP --outdir ${OUTDIR} -g 1.87e9 --bdg
