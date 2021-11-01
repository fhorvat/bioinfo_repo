#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=250g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.rnaSPAdes_ts_assembly
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=250

BASE_DIR=`basename ${PWD/rnaSPAdes*/}`
INPUT_DIR=../../../Raw/$BASE_DIR/Cleaned/trimmed
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 -name "s_*all.PE_[1,2].fq.gz"))
UNIQ_SEQ=($(printf "%s\n" "${IN_SEQ[@]%_*fq.gz}" | sort -u))
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE/_all.PE/.DN}

# ----------------Commands------------------- #
# run transcriptome assembly 
spades.py \
--rna \
-t ${THREADS} -m ${MEMORY} \
--pe-1 1 ${FILE}_1.fq.gz --pe-2 1 ${FILE}_2.fq.gz \
-o ${BASE}
