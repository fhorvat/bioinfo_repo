#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.fasta_to_2bit
#PBS -l select=ncpus=1:mem=1g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1

INPUT_DIR=.
IN_SEQ=(`find ${INPUT_DIR} -name "*.fa"`)
FILE=${IN_SEQ[0]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fa}

# script
SCRIPT=/common/WORK/fhorvat/programi/MakeHub/MakeHub/make_hub.py

# ----------------Commands------------------- #
# run
python3 ${SCRIPT} -e filip.horvat@img.cas.cz -g ${FILE} -l test_short -L test_long 
