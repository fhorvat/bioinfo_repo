#!/bin/bash
# ----------------Loading variables------------------- #
# base path and experiment 
BASE_PATH=${PWD%/Analysis/expression}
EXPERIMENT=${BASE_PATH##/*/}

# mapped and documentation
MAPPED_PATH=${BASE_PATH}/Data/Mapped/STAR_bosTau9
DOCUMENTATION_PATH=${BASE_PATH}/Data/Documentation

# genome and features
GENOME_PATH=/common/DB/genome_reference/cow/bosTau9.ARS-UCD1.2.GCA_002263795.2
FEATURES_COORDINATES=${GENOME_PATH}/ensembl.99.ARS-UCD1.2.20200415.UCSCseqnames.gtf.gz
FEATURES_NAME=${FEATURES_COORDINATES##/*/}
FEATURES_NAME=${FEATURES_NAME%.gtf.gz}
GENES_INFO_PATH=${GENOME_PATH}/${FEATURES_NAME}.geneInfo.csv

# other
SINGLE_END=TRUE
THREADS=1
GROUPING_VARIABLES="stage"
RESULTS_GROUPS="NULL"
PROTEIN_CODING_ONLY="no"
EXPLORATORY_ANALYSIS="yes"
INTERACTIVE_PLOTS="yes"
LFC_CUT="2"
PADJ_CUT="0.05"
