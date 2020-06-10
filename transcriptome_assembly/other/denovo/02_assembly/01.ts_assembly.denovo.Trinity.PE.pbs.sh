#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=250g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.trinity_ts_assembly
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=250G

BASE_DIR=`basename ${PWD/trinity/}`
INPUT_DIR=../../../Raw/$BASE_DIR/Cleaned
IN_SEQ=(`ls ${INPUT_DIR}/s_*all.PE_[1,2].txt.gz`)
UNIQ_SEQ=(`printf "%s\n" "${IN_SEQ[@]%_*txt.gz}" | sort -u`)
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}

OUT_DIR=./trinity

TRINITY_PAR="--seqType fq \
--CPU $THREADS \
--max_memory $MEMORY"

# ----------------Commands------------------- #
# run transcriptome assembly 
Trinity $TRINITY_PAR --left ${FILE}_1.txt.gz --right ${FILE}_2.txt.gz --output ${OUT_DIR} --full_cleanup

# rename .fasta from output dir 
mv ${OUT_DIR}.Trinity.fasta ${BASE}.Trinity.fasta

# get assembly statistics
/common/WORK/fhorvat/programi/Trinity/trinityrnaseq-Trinity-v2.8.4/util/TrinityStats.pl ${BASE}.Trinity.fasta > ${BASE}.Trinity.stats
