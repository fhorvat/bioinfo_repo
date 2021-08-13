#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.filter_bam
#PBS -l select=ncpus=6:mem=20g
#PBS -j oe
#PBS -J 0-4
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# number of threads
THREADS=6

# files
INPUT_DIR=/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.8.5dpp.run_2.small_RNAseq/Data/Mapped/STAR_mesAur1
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 -name "*.SE.bam"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}

BED=/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.testis/small_RNAseq/recalculate_clusters_RPM_and_RPKM/MesAur1.1k_pachytene_clusters.200730.bed

# ----------------Commands------------------- #
# get only PERFECT reads overlapping .bed
samtools view -@ ${THREADS} -b -L ${BED} ${FILE} > ./${BASE}.bam

# index
samtools index -@ ${THREADS} ${BASE}.bam
