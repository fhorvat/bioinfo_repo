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

BASE=`basename ${PWD/trinity/}`
INPUT_DIR=../../../Raw/$BASE/Cleaned

# first reads in a pair
FILES_1=(`ls ${INPUT_DIR}/*_1.txt.gz | grep -v "all.PE_1.txt.gz"`)
FILES_1=$(printf ",%s" ${FILES_1[@]})
FILES_1=${FILES_1:1}

# second reads in a pair
FILES_2=(`ls ${INPUT_DIR}/*_2.txt.gz | grep -v "all.PE_2.txt.gz"`)
FILES_2=$(printf ",%s" ${FILES_2[@]})
FILES_2=${FILES_2:1}

OUT_DIR=./trinity

TRINITY_PAR="--seqType fq \
--CPU $THREADS \
--max_memory $MEMORY"

# ----------------Commands------------------- #
# run transcriptome assembly 
Trinity $TRINITY_PAR --left $FILES_1 --right $FILES_2 --output ${OUT_DIR} --full_cleanup

# rename .fasta from output dir 
mv ${OUT_DIR}.Trinity.fasta ${BASE}.Trinity.fasta

# get assembly statistics
/common/WORK/fhorvat/programi/Trinity/trinityrnaseq-Trinity-v2.8.4/util/TrinityStats.pl ${BASE}.Trinity.fasta > ${BASE}.Trinity.stats
