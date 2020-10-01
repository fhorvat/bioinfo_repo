#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=50g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.00.unpigz_fastq
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=50G

IN_DIR=..
IN_SEQ=(${IN_DIR}/s_mouse*.PE_1.txt.gz)
UNIQ_SEQ=($(printf "%s\n" "${IN_SEQ[@]%_*.txt.gz}" | sort -u))
FILE=${UNIQ_SEQ[0]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.txt.gz}

# ----------------Commands------------------- #
# unpigz
unpigz -p $THREADS -c ${FILE}_1.txt.gz > ${BASE}_1.fastq
unpigz -p $THREADS -c ${FILE}_2.txt.gz > ${BASE}_2.fastq
