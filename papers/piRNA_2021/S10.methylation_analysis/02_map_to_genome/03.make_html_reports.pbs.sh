#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03.make_html_reports
#PBS -l select=ncpus=1:mem=5g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# source
source ./000.load_variables.sh

# override number of threads
THREADS=1
MEMORY=5G

# input
INPUT_DIR=.
IN_SEQ=($(find . ${INPUT_DIR} -maxdepth 1 -name "*_bismark_bt2_SE_report.txt"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%_bismark_bt2_SE_report.txt}

# ----------------Commands------------------- #
# run 
bismark2report
bismark2summary
