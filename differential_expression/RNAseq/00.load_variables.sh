#!/bin/bash
# ----------------Loading variables------------------- #
# variables to change - mapped dir and genome
MAPPED_DIR="STAR_mm10"
IN_GENOME=${MOUSE}
ENSEMBL_VERSION=99

# more variables to change 
GROUPING_VARIABLES="genotype"
RESULTS_GROUPS="Mov10l_KO,Mov10l_WT Mov10l_HET,Mov10l_WT Mov10l_KO,Mov10l_HET" # "no" for no diff. exp. analysis
SINGLE_END=TRUE

# even more variables to change, but usually left alone 
PROTEIN_CODING_ONLY="no"
EXPLORATORY_ANALYSIS="yes"
VULCANO_PLOTS="no"
INTERACTIVE_PLOTS="no"
LFC_CUT="0.5"
PADJ_CUT="0.05"
FPKM_CUT="0"
THREADS=1

# base path and experiment 
BASE_PATH=${PWD%/Analysis/*}
EXPERIMENT=${BASE_PATH##/*/}

# mapped and documentation
MAPPED_PATH=${BASE_PATH}/Data/Mapped/${MAPPED_DIR}
DOCUMENTATION_PATH=${BASE_PATH}/Data/Documentation

# genome and features
GENOME_BASE=/common/DB/genome_reference
GENOME_PATH=${GENOME_BASE}/${IN_GENOME}
FEATURES_COORDINATES=$(find ${GENOME_PATH} -maxdepth 1 -name "ensembl.${ENSEMBL_VERSION}.*.UCSCseqnames.gtf.gz")
FEATURES_NAME=${FEATURES_COORDINATES##/*/}
FEATURES_NAME=${FEATURES_NAME%.gtf.gz}
GENES_INFO_PATH=${GENOME_PATH}/${FEATURES_NAME}.geneInfo.csv
