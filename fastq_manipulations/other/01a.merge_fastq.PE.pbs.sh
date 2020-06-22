#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=30g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01a.merge_fastq.PE
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# input dir
IN_DIR=.

# base name
FILES_1=($(find ${IN_DIR} -name "*_1.txt.gz"))
BASE=${FILES_1[0]}
BASE=${BASE%_r*.txt.gz}
BASE=${BASE#${IN_DIR}/}

# first reads in pair
FILES_1=$(echo ${FILES_1[@]} | paste -sd" " -) 

# second reads in a pair
FILES_2=$(find ${IN_DIR} -name "${BASE}*_2.txt.gz" | paste -sd" " -)

# ----------------Commands------------------- #
# concatenate
cat ${FILES_1} > ${BASE}_all.PE_1.txt.gz
cat ${FILES_2} > ${BASE}_all.PE_2.txt.gz
