#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=2:mem=10g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N bbduk
#PBS -J 0-15
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
ADAPTER=/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/common/Yang_2016_development_smallRNA_PRJNA257532/Data/documentation/adapter.fa 

IN_DIR=/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/common/Yang_2016_development_smallRNA_PRJNA257532/Data/Raw/Links
IN_SEQ=($IN_DIR/*.fastq.gz)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE##*$IN_DIR/}
BASE=${BASE%%.fastq.gz}

BBDUK_PAR="overwrite=t \
k=21 \
mink=12 \
hdist=1 \
threads=1 \
qtrim=r \
minlength=15 \
minoverlap=8"

# ----------------Commands------------------- #
bbduk2.sh in=$FILE rref=$ADAPTER out=${BASE}_cl.fastq.gz stats=${BASE}_cl.log $BBDUK_PAR 2> ${BASE}_cl.stats
