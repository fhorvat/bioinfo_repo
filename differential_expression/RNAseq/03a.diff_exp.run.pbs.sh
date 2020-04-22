#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03.diff_exp
#PBS -l select=ncpus=1:mem=20g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script path
SCRIPT=03b.diff_exp.script.R

# read variables from source
source ./00.load_variables.sh

# counts path
COUNTS_PATH=($(find . -maxdepth 1 -name "*counts.txt"))

# input for manual script
echo -e '\n'\
experiment=\'$EXPERIMENT\''\n'\
single_end=$SINGLE_END'\n'\
threads=$THREADS'\n'\
mapped_path=\'$MAPPED_PATH\''\n'\
documentation_path=\'$DOCUMENTATION_PATH\''\n'\
features_coordinates=\'$FEATURES_COORDINATES\''\n'\
features_name=\'$FEATURES_NAME\''\n'\
genes_info_path=\'$GENES_INFO_PATH\''\n'\
grouping_variables=\'$GROUPING_VARIABLES\''\n'\
results_groups=\'$RESULTS_GROUPS\''\n'\
protein_coding_only=\'$PROTEIN_CODING_ONLY\''\n'\
exploratory_analysis=\'$EXPLORATORY_ANALYSIS\''\n'\
interactive_plots=\'$INTERACTIVE_PLOTS\''\n'\
counts_path=\'$COUNTS_PATH\''\n'\
lfc_cut=\'$LFC_CUT\''\n'\
padj_cut=\'$PADJ_CUT\''\n'\
'\n'\
'results_groups <- str_split(results_groups, pattern = " ") %>% unlist()''\n'\
'grouping_variables <- str_split(grouping_variables, pattern = " ") %>% unlist()''\n'

# ----------------Commands------------------- #
# run script
Rscript $SCRIPT \
--experiment $EXPERIMENT \
--single_end $SINGLE_END \
--threads $THREADS \
--mapped_path $MAPPED_PATH \
--documentation_path $DOCUMENTATION_PATH \
--features_coordinates $FEATURES_COORDINATES \
--features_name $FEATURES_NAME \
--genes_info_path $GENES_INFO_PATH \
--grouping_variables $GROUPING_VARIABLES \
--results_groups $RESULTS_GROUPS \
--protein_coding_only $PROTEIN_CODING_ONLY \
--exploratory_analysis $EXPLORATORY_ANALYSIS \
--interactive_plots $INTERACTIVE_PLOTS \
--counts_path $COUNTS_PATH \
--lfc_cut $LFC_CUT \
--padj_cut $PADJ_CUT
