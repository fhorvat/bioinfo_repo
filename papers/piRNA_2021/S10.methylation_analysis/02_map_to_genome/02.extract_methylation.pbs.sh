#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.extract_methylation
#PBS -l select=ncpus=10:mem=20g
#PBS -j oe
#PBS -J 0-2
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# source
source ./000.load_variables.sh

# override number of threads
THREADS=3
MEMORY=20G

# input
IN_SEQ=($(find . -maxdepth 1 -name "*.deduplicated.bam"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}

# script parameters
SCRIPT_PAR="--parallel ${THREADS} \
--comprehensive \
--bedGraph \
--gzip \
--buffer_size ${MEMORY} \
--output ."

# ----------------Commands------------------- #
# run 
bismark_methylation_extractor ${SCRIPT_PAR} ${FILE}
