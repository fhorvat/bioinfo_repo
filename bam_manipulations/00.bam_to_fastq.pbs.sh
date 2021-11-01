#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.00.bam_to_fastq
#PBS -J 0-3
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12

IN_DIR=/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.smallRNAseq/Data/Mapped/STAR_Siomi.multimappers/6_filter_18to32nt
IN_SEQ=($(find $IN_DIR -maxdepth 1 \( -name "*.bam" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.bam}

# ----------------Commands------------------- #
# get mapped reads to fastq
samtools fastq -@ $THREADS ${FILE} > ${BASE}.fastq

