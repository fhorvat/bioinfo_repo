#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.stranded_coverage_clusters
#PBS -l select=ncpus=1:mem=10g
#PBS -j oe
#PBS -J 0-5
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script path
SCRIPT=01b.stranded_coverage_clusters.script.R

# clusters table path
INPUT_DIR=..
CLUSTERS_PATH=$(find $INPUT_DIR -name "MesAur1.1k_pachytene_clusters.200730.recalculated.fused.PS.xlsx")

# bam files path
INPUT_DIR=/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq.reseq/Data/Mapped/STAR_mesAur1
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 \( -name "*.SE.24to31nt.bam" -not -name "*HET*" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}

# input for manual script
echo -e '\n'\
clusters_path=\'$CLUSTERS_PATH\''\n'\
bam_path=\'$FILE\''\n'

# ----------------Commands------------------- #
# run script
Rscript $SCRIPT \
--clusters_path $CLUSTERS_PATH \
--bam_path $FILE
