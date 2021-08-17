#!/bin/bash
# ----------------Loading variables------------------- #
# set experiment
EXPERIMENT=$(basename $PWD)
EXPERIMENT_NAME=${EXPERIMENT}".miRNA"

# mapped and documentation
MAPPED_PATH=../../datasets/${EXPERIMENT}/5_perfect_reads
DOCUMENTATION_PATH=../../Documentation/${EXPERIMENT}

# genome and features
FEATURES_COORDINATES=../../annotation/${EXPERIMENT}
FEATURES_NAME=${EXPERIMENT_NAME}
GENES_INFO_PATH=${FEATURES_COORDINATES}.geneInfo.csv

# other
THREADS=1
GROUPING_VARIABLES="stage"
RESULTS_GROUPS="no"
EXPLORATORY_ANALYSIS="yes"
INTERACTIVE_PLOTS="no"
