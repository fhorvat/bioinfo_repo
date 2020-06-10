#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=250g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.ts_assembly.trinity
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=250G

INPUT_DIR=..
IN_SEQ=(`ls ${INPUT_DIR}/*.txt.gz`)
FILES=$(printf ",%s" ${IN_SEQ[@]})
FILES=${FILES:1}

BASE="pwd"
OUT_DIR=${BASE}_trinity_assembly

TRINITY_PAR="--seqType fq \
--CPU $THREADS \
--max_memory $MEMORY"

# ----------------Commands------------------- #
# run transcriptome assembly 
Trinity $TRINITY_PAR --single $FILES --output ${OUT_DIR} --full_cleanup

# rename .fasta from output dir 
mv ${OUT_DIR}.Trinity.fasta ${BASE}.Trinity.fasta

# get assembly statistics
/common/WORK/fhorvat/programi/Trinity/trinityrnaseq-Trinity-v2.5.1/util/TrinityStats.pl ${BASE}.Trinity.fasta > ${BASE}.Trinity.stats
