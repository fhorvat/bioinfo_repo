#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=250g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03.oases_ts_assembly
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=250G
OMP_THREAD_LIMIT=$THREADS

DIR_BASE=`basename ${PWD/oases/}`
INPUT_DIR=../../../Raw/$DIR_BASE/Cleaned
IN_SEQ=(`ls ${INPUT_DIR}/*all.PE_[1,2].txt.gz`)
UNIQ_SEQ=(`printf "%s\n" "${IN_SEQ[@]%_[1-2].txt.gz}" | sort -u`)
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.txt.gz}

OASES_SCRIPT=/common/WORK/fhorvat/programi/velvet_oases/oases/scripts/oases_pipeline.py

# ----------------Commands------------------- #
# get insert size
HIST_FILE=(`ls *.hist.txt`)
INSERT_SIZE=(`sed -n 1p $HIST_FILE | awk '{print $2}' | xargs printf "%.*f\n" 0`)

# run transcriptome assembly 
$OASES_SCRIPT \
-d " -fastq.gz -shortPaired ${FILE}_1.txt.gz ${FILE}_2.txt.gz " \
-m 21 \
-M 35 \
-o $OUTPUT_DIR \
-p " -ins_length $INSERT_SIZE -min_trans_lgth 100 "\
-c \
--merge=27

