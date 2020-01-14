#!/bin/bash
# ----------------Loading variables------------------- #
# set experiment
EXPERIMENT=DicerX_embryos

# base path
BASE_PATH=/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/DicerX_viral_infection/datasets/2019_Dec

# mapped and documentation
MAPPED_PATH=${BASE_PATH}/Data/Mapped/STAR_mm10
DOCUMENTATION_PATH=${BASE_PATH}/Data/Documentation

# genome and features
GENOME_PATH=/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2
FEATURES_COORDINATES=${GENOME_PATH}/miRBase.22.mm10.20181605.gff3
FEATURES_NAME=${FEATURES_COORDINATES##/*/}
FEATURES_NAME=${FEATURES_NAME%.gff3}

# other
THREADS=1
GROUPING_VARIABLES="genotype"
RESULTS_GROUPS="KO,WT HET,WT KO,HET"
EXPLORATORY_ANALYSIS="yes"
INTERACTIVE_PLOTS="yes"
