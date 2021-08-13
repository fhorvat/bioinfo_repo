#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.06.coverage2cytosine
#PBS -l select=ncpus=1:mem=10g
#PBS -J 0-2
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# source
source ./000.load_variables.sh

# override number of threads
THREADS=1

# input 
INPUT_DIR=.
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 \( -name "*.bismark.cov" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bismark.cov}

# genome folder
GENOME_PATH=../../../bismark_index/FLI_consensus.IAPLTR3

# ----------------Commands------------------- #
# get cytosine methylation
coverage2cytosine --merge_CpG --genome_folder $GENOME_PATH -o ${BASE} ${FILE}

