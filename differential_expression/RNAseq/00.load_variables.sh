#!/bin/bash
# ----------------Loading variables------------------- #
# base path and experiment 
BASE_PATH=${PWD%/Analysis/*}
EXPERIMENT=${BASE_PATH##/*/}

# mapped and documentation
MAPPED_PATH=${BASE_PATH}/Data/Mapped/STAR_mm10
DOCUMENTATION_PATH=${BASE_PATH}/Data/Documentation

# genome and features
GENOME_PATH=/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2
FEATURES_COORDINATES=${GENOME_PATH}/ensembl.99.GRCm38.p6.20200415.UCSCseqnames.gtf.gz
FEATURES_NAME=${FEATURES_COORDINATES##/*/}
FEATURES_NAME=${FEATURES_NAME%.gtf.gz}
GENES_INFO_PATH=${GENOME_PATH}/${FEATURES_NAME}.geneInfo.csv

# other
SINGLE_END=TRUE
THREADS=1
GROUPING_VARIABLES="genotype"
RESULTS_GROUPS="Mov10l_KO,Mov10l_WT Mov10l_HET,Mov10l_WT Mov10l_KO,Mov10l_HET" # "no" for no diff. exp. analysis 
PROTEIN_CODING_ONLY="no"
EXPLORATORY_ANALYSIS="yes"
VULCANO_PLOTS="no"
INTERACTIVE_PLOTS="no"
LFC_CUT="0.5"
PADJ_CUT="0.1"
FPKM_CUT="1"
