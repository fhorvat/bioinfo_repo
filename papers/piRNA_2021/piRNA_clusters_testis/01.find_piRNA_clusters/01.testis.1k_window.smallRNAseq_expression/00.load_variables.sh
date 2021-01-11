#!/bin/bash
# ----------------Loading variables------------------- #
# set experiment
EXPERIMENT=Mov10l1_KO

# base path
BASE_PATH=/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq

# mapped and documentation
MAPPED_PATH=${BASE_PATH}/Data/Mapped/STAR_mesAur1.new
DOCUMENTATION_PATH=${BASE_PATH}/Data/Documentation

# genome and features
GENOME_PATH=/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1
FEATURES_COORDINATES=MesAur1.1k_windows.gff
FEATURES_NAME=${FEATURES_COORDINATES##/*/}
FEATURES_NAME=${FEATURES_NAME%.gff}

# other
SINGLE_END=TRUE
THREADS=1
GROUPING_VARIABLES="genotype age"
RESULTS_GROUPS="Mov10l_KO_13dpp,Mov10l_WT_13dpp"
PROTEIN_CODING_ONLY="no"
EXPLORATORY_ANALYSIS="yes"
INTERACTIVE_PLOTS="no"
