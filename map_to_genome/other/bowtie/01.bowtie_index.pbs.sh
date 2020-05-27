#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=50g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.bowtie_index
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script parameters
THREADS=12

# set input variables
IN_DIR=.
IN_GENOME=($(find $IN_DIR -maxdepth 1 \( -name "*fa" -not -name "ensembl*" \)))
GENOME_NAME=`basename $IN_GENOME`
GENOME_NAME=${GENOME_NAME%.fa}.bowtie1_index

# ----------------Commands------------------- #
# generate index
/common/WORK/fhorvat/programi/bowtie/bowtie-1.2.1.1/bowtie-build --threads $THREADS $IN_GENOME $GENOME_NAME 
