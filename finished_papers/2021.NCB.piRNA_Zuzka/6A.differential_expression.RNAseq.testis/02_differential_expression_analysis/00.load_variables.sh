#!/bin/bash
# ----------------Loading variables------------------- #
# base path and experiment 
BASE_PATH=${PWD%/Analysis/*}
EXPERIMENT=${BASE_PATH##/*/}

# mapped and documentation
MAPPED_PATH=${BASE_PATH}/Data/Mapped/STAR_mesAur1
DOCUMENTATION_PATH=${BASE_PATH}/Data/Documentation

# genome and features
GENOME_PATH=/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1/gtf_with_some_refSeq_genes
FEATURES_COORDINATES=${GENOME_PATH}/ensembl.99.MesAur1.0.20200415.UCSCseqnames.Piwil3_Tex101.from_RefSeq.gtf.gz
FEATURES_NAME=${FEATURES_COORDINATES##/*/}
FEATURES_NAME=${FEATURES_NAME%.gtf.gz}
GENES_INFO_PATH=${GENOME_PATH}/${FEATURES_NAME}.geneInfo.csv

# other
SINGLE_END=TRUE
THREADS=1
GROUPING_VARIABLES="genotype age"
RESULTS_GROUPS="Mov10l1_KO_8.5dpp,Mov10l1_WT_8.5dpp" # "no" for no diff. exp. analysis 
PROTEIN_CODING_ONLY="yes"
EXPLORATORY_ANALYSIS="yes"
VULCANO_PLOTS="no"
INTERACTIVE_PLOTS="yes"
LFC_CUT="0.5"
PADJ_CUT="0.01"
FPKM_CUT="0"
