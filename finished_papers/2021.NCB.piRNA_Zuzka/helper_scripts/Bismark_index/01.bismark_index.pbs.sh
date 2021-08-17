#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=24:mem=50g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.bismark_index
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script parameters
THREADS=12

# set input variables
INPUT_DIR=.
#IN_SEQ=$(find $INPUT_DIR -name "mm10.fa.gz" -or -name "mm10.fa")
#FILE=${IN_SEQ[0]}
#BASE=${BASE#${INPUT_DIR}/}
#BASE=${BASE%.fa.gz}
#BASE=${BASE%.fa}

# ----------------Commands------------------- #
# generate index
bismark_genome_preparation --verbose --bowtie2 --parallel ${THREADS} ${INPUT_DIR}
