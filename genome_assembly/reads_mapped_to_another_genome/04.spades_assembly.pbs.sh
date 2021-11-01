#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=40g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.04.spades_assembly
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=40

IN_DIR=.

IN_FRAGMENT=($(find ${IN_DIR} -maxdepth 1 \( -name "*_s.fastq" -and -name "*fragment*"  \)))
NAME_FRAGMENT=${IN_FRAGMENT%_s.fastq}

IN_MATE=($(find ${IN_DIR} -maxdepth 1 \( -name "*_s.fastq" -and -name "*3kb*"  \)))
NAME_MATE=${IN_MATE%_s.fastq}

BASE=${NAME_FRAGMENT/_fragment/}.assembly

# ----------------Commands------------------- #
# assembly
spades.py \
-o ${BASE} \
--only-assembler \
-t ${THREADS} -m ${MEMORY} \
-k 33,55,77 \
--pe-1 1 ${NAME_FRAGMENT}_1.fastq --pe-2 1 ${NAME_FRAGMENT}_2.fastq --pe-s 1 ${NAME_FRAGMENT}_s.fastq \
--mp-1 1 ${NAME_MATE}_1.fastq --mp-2 1 ${NAME_MATE}_2.fastq --mp-s 1 ${NAME_MATE}_s.fastq 
