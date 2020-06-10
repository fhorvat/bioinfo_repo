#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=100g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.ts_assembly.trinity
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=100G

INPUT_DIR=.
FILE=(`ls ${INPUT_DIR}/*Aligned.sortedByCoord.out.bam`)
BASE="pwd"
OUT_DIR=${BASE}_trinity_assembly

TRINITY_PAR="--genome_guided_max_intron 1000000 \
--full_cleanup \
--CPU $THREADS \
--max_memory $MEMORY"

# ----------------Commands------------------- #
# run transcriptome assembly 
Trinity $TRINITY_PAR --genome_guided_bam $FILE --output ${OUT_DIR}

# move .fasta from output dir 
mv ${OUT_DIR}/Trinity-GG.fasta .
mv Trinity-GG.fasta ${BASE}.Trinity.fasta

# remove output dir
[ -f "${BASE}.Trinity.fasta" ] && rm -r ${OUT_DIR}

# get assembly statistics
/common/WORK/fhorvat/programi/Trinity/trinityrnaseq-Trinity-v2.5.1/util/TrinityStats.pl ${BASE}.Trinity.fasta > ${BASE}.Trinity.stats
