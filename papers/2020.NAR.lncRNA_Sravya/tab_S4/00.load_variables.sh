#!/bin/bash
# ----------------Loading variables------------------- #
# set experiment
EXPERIMENT=Fugaku.perfect

# base path
BASE_PATH=/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Fugaku/Data/Mapped/STAR_mm10.perfect

# mapped and documentation
MAPPED_PATH=/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Fugaku/Data/Mapped/STAR_mm10.perfect
DOCUMENTATION_PATH=/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Fugaku/Data/Documentation

# genome and features
GENOME_PATH=/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2
FEATURES_COORDINATES=/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/pseudogenes_expression/2019_Ganesh_Table_S4_Elob_pseudogenes.saf
FEATURES_NAME=${FEATURES_COORDINATES##/*/}
FEATURES_NAME=${FEATURES_NAME%.saf}
GENES_INFO_PATH=${GENOME_PATH}/ensembl.93.GRCm38.p6.20180919.UCSCseqnames.geneInfo.csv

# other
SINGLE_END=FALSE
THREADS=1
GROUPING_VARIABLES="stage"
RESULTS_GROUPS="GV.WE"
PROTEIN_CODING_ONLY="no"
EXPLORATORY_ANALYSIS="yes"
INTERACTIVE_PLOTS="yes"
