#!/bin/bash
# ----------------Loading variables------------------- #
# set experiment
EXPERIMENT=hamster_oocyte_Mov10l.smallRNAseq

# base path
BASE_PATH=/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.deduplicated.smallRNAseq

# mapped and documentation
MAPPED_PATH=../split_by_read_length
DOCUMENTATION_PATH=${BASE_PATH}/Data/Documentation

# genome and features
GENOME_PATH=/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1
FEATURES_COORDINATES=$(find ../ -name "piRNA_clusters.*gff")
FEATURES_NAME=${FEATURES_COORDINATES##*/}
FEATURES_NAME=${FEATURES_NAME%.gff}

# other
SINGLE_END=TRUE
THREADS=1
GROUPING_VARIABLES="genotype"
RESULTS_GROUPS="Mov10l1_WT Mov10l1_HET Mov10l1_KO"
PROTEIN_CODING_ONLY="no"
EXPLORATORY_ANALYSIS="yes"
INTERACTIVE_PLOTS="no"
