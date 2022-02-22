#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02a.plot_repeat_landscape
#PBS -l select=ncpus=1:mem=1g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script path
SCRIPT=02b.plot_repeat_landscape.script.R

# divergence file
INPUT_DIR=.
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 -name "*.divsum"))
FILE=${IN_SEQ[0]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.divsum}

# 2bit genome file
INPUT_DIR=/common/DB/genome_reference/${SIOMI}
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 -name "*.2bit"))
FILE_GENOME=${IN_SEQ[0]}

# hardcoded repeat classes and colors
SCRIPT_PATH=/common/WORK/fhorvat/Projekti/Svoboda/scripts/repeatMasker/repeat_landscape
CLASSES_PATH=${SCRIPT_PATH}/rmsk_landscape.class_names_to_graph_names.txt
COLORS_PATH=${SCRIPT_PATH}/rmsk_landscape.graph_names_colors.txt

# input for manual script
echo -e '\n'\
div_path=\'${FILE}\''\n'\
genome_2bit_path=\'${FILE_GENOME}\''\n'\
classes_path=\'${CLASSES_PATH}\''\n'\
colors_path=\'${COLORS_PATH}\''\n'\
plot_name=\'${BASE}.repeat_landscape\''\n'

# ----------------Commands------------------- #
# run script
Rscript ${SCRIPT} \
--div_path ${FILE} \
--genome_2bit_path ${FILE_GENOME} \
--classes_path ${CLASSES_PATH} \
--colors_path ${COLORS_PATH} \
--plot_name ${BASE}.repeat_landscape
