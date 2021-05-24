#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.00.join_genomes
#PBS -l select=ncpus=6:mem=10g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# threads and memory
THREADS=6

# input and output
INPUT_DIR=..
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 -type f -name "mm10.fa"))
FILE_1=${IN_SEQ[0]}
BASE_1=${FILE_1#${INPUT_DIR}/}
BASE_1=${BASE_1%.fa}

INPUT_DIR=/common/WORK/fhorvat/Projekti/Svoboda/antiviral_RNAi.Marcos/datasets/Documentation
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 -type f -name "*.fasta"))
FILE_2=$(printf "%s " ${IN_SEQ[@]})
BASE_2="TBEV.LMCV"

# ----------------Commands------------------- #
# join together
cat ${FILE_1} ${FILE_2} > ${BASE_1}.${BASE_2}.fa

# zip
bgzip -@ $THREADS ${BASE_1}.${BASE_2}.fa

# index
samtools faidx ${BASE_1}.${BASE_2}.fa.gz
