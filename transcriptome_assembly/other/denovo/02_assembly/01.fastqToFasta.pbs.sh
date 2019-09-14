#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=16:mem=30g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.fastq_to_fasta
#PBS -J 0-1
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=16
DIR_BASE=`basename ${PWD/oases/}`

INPUT_DIR=../../../Raw/$DIR_BASE/Cleaned
IN_SEQ=(`ls ${INPUT_DIR}/*all.PE_[1,2].txt.gz`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.txt.gz}

# ----------------Commands------------------- #
# trasform fastq to fasta
unpigz -cp $THREADS $FILE | sed -n '1~4s/^@/>/p;2~4p' > ${BASE}.fasta
