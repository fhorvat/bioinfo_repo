#!/bin/bash
# ----------------Loading variables------------------- #
# set experiment
EXPERIMENT=hamster_Siomi_WT.Piwil1_Piwil3_IP.small_RNAseq

# base path
BASE_PATH=/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_Siomi_WT.Piwil1_Piwil3_IP.small_RNAseq

# mapped and documentation
MAPPED_PATH=${BASE_PATH}/Data/Mapped/STAR_mesAur1
DOCUMENTATION_PATH=${BASE_PATH}/Data/Documentation

# genome and features
GENOME_PATH=/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1
FEATURES_COORDINATES=$(find ./ -name "piRNA_clusters.oocyte_PIWIL1.rpkm_cutoff.1.20210513.gff")
FEATURES_NAME=${FEATURES_COORDINATES##*/}
FEATURES_NAME=${FEATURES_NAME%.gff}

# other
SINGLE_END=TRUE
THREADS=1
GROUPING_VARIABLES="IP"
RESULTS_GROUPS="Mov10l_KO_13dpp,Mov10l_WT_13dpp"
PROTEIN_CODING_ONLY="no"
EXPLORATORY_ANALYSIS="yes"
INTERACTIVE_PLOTS="no"
