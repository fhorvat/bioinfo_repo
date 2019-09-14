#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=30:mem=250g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.bosTau.trinity.gg
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=30
MEMORY=250G

BASE_DIR=`basename ${PWD/trinity*/}`
INPUT_DIR=../../../Mapped/$BASE_DIR/4_merged_replicates
IN_SEQ=(`ls ${INPUT_DIR}/s_*all.PE.bam`)
FILE=${IN_SEQ[0]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE/_all.PE.bam/.GG}

OUT_DIR=./trinity

TRINITY_PAR="--genome_guided_max_intron 589824 \
--CPU $THREADS \
--max_memory $MEMORY"

# ----------------Commands------------------- #
# run transcriptome assembly 
Trinity $TRINITY_PAR --genome_guided_bam ${FILE} --output ${OUT_DIR} --full_cleanup

# move .fasta from output dir 
mv ${OUT_DIR}/Trinity-GG.fasta ./${BASE}.Trinity.fasta

# get assembly statistics
/common/WORK/fhorvat/programi/Trinity/trinityrnaseq-Trinity-v2.8.4/util/TrinityStats.pl ${BASE}.Trinity.fasta > ${BASE}.Trinity.fasta.stats
