#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=2:mem=10g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N skewer
#PBS -J 0-15
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
ADAPTER=/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/common/Yang_2016_development_smallRNA_PRJNA257532/Data/documentation/adapter_revcomp.fa
TRIMMOMATIC='/common/WORK/fhorvat/programi/Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar'
THREADS=2

IN_DIR=/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/common/Yang_2016_development_smallRNA_PRJNA257532/Data/Raw/Links
IN_SEQ=($IN_DIR/*.fastq.gz)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE##*$IN_DIR/}
BASE=${BASE%%.fastq.gz}

TRIM_PAR="-l 15
-z
-t $THREADS"

# ----------------Commands------------------- #
skewer $FILE -x $ADAPTER -o ${BASE} $TRIM_PAR
