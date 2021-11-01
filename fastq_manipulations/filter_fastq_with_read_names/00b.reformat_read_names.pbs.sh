#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.00b.reformat_read_names
#PBS -j oe
#PBS -J 0-3
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1
MEMORY=20G

# input fastq file
IN_DIR=..
IN_SEQ=($(find ${IN_DIR} \( -name "*.atrim.txt.gz" -and -name "s*oocytes*" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.atrim.txt.gz}

# ----------------Commands------------------- #
# remove everything after space in read names
zcat $FILE | sed -e 's/ .*//g' > ${BASE}.fastq
