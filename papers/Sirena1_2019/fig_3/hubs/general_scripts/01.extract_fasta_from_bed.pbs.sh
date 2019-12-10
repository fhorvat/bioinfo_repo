#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.fasta_from_bed
#PBS -l select=ncpus=1:mem=1g
#PBS -J 0-10
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1

INPUT_DIR=../coordinates
IN_SEQ=(`find ${INPUT_DIR} -name "*.inter_exons.bed"`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bed}
SCAFFOLD=`awk '{print $1}' $FILE`

INPUT_DIR_FASTA=/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/rodents_evolution/genomes/unpacked
GENOME_FASTA=($INPUT_DIR_FASTA/${BASE%.inter_exons}*.fa)

# ----------------Commands------------------- #
# get fasta on correct strand with bedtools
samtools faidx $GENOME_FASTA $SCAFFOLD > ${BASE%.inter_exons}.${SCAFFOLD}.fa
