#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03.diff_exp
#PBS -l select=ncpus=1:mem=40g
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
threads=$THREADS'\n'\
mapped_path=\'$MAPPED_PATH\''\n'\
documentation_path=\'$DOCUMENTATION_PATH\''\n'\
features_name=\'$FEATURES_NAME\''\n'\
grouping_variables=\'$GROUPING_VARIABLES\''\n'\
results_groups=\'$RESULTS_GROUPS\''\n'\
exploratory_analysis=\'$EXPLORATORY_ANALYSIS\''\n'\
interactive_plots=\'$INTERACTIVE_PLOTS\''\n'\
counts_path=\'$COUNTS_PATH\''\n'\
'\n'\
'results_groups <- str_split(results_groups, pattern = " ") %>% unlist()''\n'\
'grouping_variables <- str_split(grouping_variables, pattern = " ") %>% unlist()''\n'

# ----------------Commands------------------- #
# run script
Rscript $SCRIPT \
--experiment $EXPERIMENT \
--threads $THREADS \
--mapped_path $MAPPED_PATH \
--documentation_path $DOCUMENTATION_PATH \
--features_name $FEATURES_NAME \
--grouping_variables $GROUPING_VARIABLES \
--results_groups $RESULTS_GROUPS \
--exploratory_analysis $EXPLORATORY_ANALYSIS \
--interactive_plots $INTERACTIVE_PLOTS \
--counts_path $COUNTS_PATH
