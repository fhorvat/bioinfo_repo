#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=30g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.04.merge_fastq.PE
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
IN_DIR=.
# first reads in a pair
FILES_1=(`ls ${IN_DIR}/*_1.clean.trim.txt.gz`)
BASE=${FILES_1[0]}
BASE=${BASE%_r*.clean.trim.txt.gz}
BASE=${BASE#${IN_DIR}/}

# second reads in a pair
FILES_2=(`ls ${IN_DIR}/*_2.clean.trim.txt.gz`)

# ----------------Commands------------------- #
# concatenate
cat `echo ${FILES_1[@]}` > ${BASE}_all.PE_1.txt.gz
cat `echo ${FILES_2[@]}` > ${BASE}_all.PE_2.txt.gz
