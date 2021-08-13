#!/bin/bash
# ----------------Loading variables------------------- #
# set experiment
EXPERIMENT=Mov10l1_KO

# base path
BASE_PATH=/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.8.5dpp.small_RNAseq

# mapped and documentation
MAPPED_PATH=${BASE_PATH}/Data/Mapped/STAR_mesAur1
DOCUMENTATION_PATH=${BASE_PATH}/Data/Documentation

# genome and features
GENOME_PATH=/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.testis/small_RNAseq/recalculate_clusters_RPM_and_RPKM
FEATURES_COORDINATES=${GENOME_PATH}/MesAur1.1k_pachytene_clusters.200730.gff
FEATURES_NAME=${FEATURES_COORDINATES##/*/}
FEATURES_NAME=${FEATURES_NAME%.gff}

# other
SINGLE_END=TRUE
THREADS=1
GROUPING_VARIABLES="genotype age"
RESULTS_GROUPS="Mov10l_KO_8.5,Mov10l_WT_8.5"
PROTEIN_CODING_ONLY="no"
EXPLORATORY_ANALYSIS="yes"
INTERACTIVE_PLOTS="no"
