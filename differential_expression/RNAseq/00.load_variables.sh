#!/bin/bash
# ----------------Loading variables------------------- #
# base path and experiment 
BASE_PATH=${PWD%/Analysis/expression}
EXPERIMENT=${BASE_PATH##/*/}

# mapped and documentation
MAPPED_PATH=${BASE_PATH}/Data/Mapped/STAR_mm10
DOCUMENTATION_PATH=${BASE_PATH}/Data/Documentation

# genome and features
GENOME_PATH=/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2
FEATURES_COORDINATES=${GENOME_PATH}/ensembl.93.GRCm38.p6.20180919.UCSCseqnames.gtf.gz
FEATURES_NAME=${FEATURES_COORDINATES##/*/}
FEATURES_NAME=${FEATURES_NAME%.gtf.gz}
GENES_INFO_PATH=${GENOME_PATH}/${FEATURES_NAME}.geneInfo.csv

# other
SINGLE_END=TRUE
THREADS=1
GROUPING_VARIABLES="genotype"
RESULTS_GROUPS="Ago2_KO,Ago2_WT Ago2_HET,Ago2_WT Ago2_KO,Ago2_HET"
PROTEIN_CODING_ONLY="no"
EXPLORATORY_ANALYSIS="yes"
INTERACTIVE_PLOTS="yes"
LFC_CUT="2"
PADJ_CUT="0.05"
