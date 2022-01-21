#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.repeat_landscape
#PBS -l select=ncpus=1:mem=1g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# input align file
INPUT_DIR=..
IN_SEQ=($(find ${INPUT_DIR} -name "*align"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.align}

# input genome 2bit
INPUT_DIR=../..
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 -name "*.2bit"))
FILE_GENOME=${IN_SEQ[$PBS_ARRAY_INDEX]}

# input scripts
SCRIPT_1=/common/WORK/fhorvat/programi/RepeatMasker/RepeatMasker-open-4-0-9-p2/util/calcDivergenceFromAlign.pl
SCRIPT_2=/common/WORK/fhorvat/programi/RepeatMasker/RepeatMasker-open-4-0-9-p2/util/createRepeatLandscape.pl

# ----------------Commands------------------- #
# create divergence from repeatMasker .align file
perl ${SCRIPT_1} -s ${BASE}.divsum ${FILE}\

# plot repeat landscape in a .html file
perl ${SCRIPT_2} -div ${BASE}.divsum -twoBit ${FILE_GENOME} > ${BASE}.repeat_landscape.html
