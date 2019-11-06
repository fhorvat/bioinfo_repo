#!/bin/bash
# ----------------Loading variables------------------- #
# get experiment
L1_LIST=${PWD##/*/}
IN_PATH=${PWD%/expression/${L1_LIST}}
EXPERIMENT=${IN_PATH##/*/}

# is single end
SINGLE_END=FALSE

# dataset path
BASE_PATH=/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression
DATASET_PATH=${BASE_PATH}/datasets/${EXPERIMENT}/Mapped/perfect_alignments.all_multimappers

# documentation path
DOCUMENTATION_PATH=/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/Documentation/${L1_LIST}
