#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.merge_bams
#PBS -l select=ncpus=12:mem=30g
#PBS -J 0-%N_JOBS
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# source
source ./000.load_variables.sh

# override number of threads
THREADS=12

# input
REPLICATE=${REPLICATES[$PBS_ARRAY_INDEX]}
FILES_LIST=($(find ${INPUT_DIR} -maxdepth 1 \( -name "${REPLICATE}*r*.*E.*bam" \)))
FILES_LIST=$(printf "%s " ${FILES_LIST[@]})

# just in case
echo ${FILES_LIST}

# ----------------Commands------------------- #
# merge
samtools merge ${REPLICATE}.${LAYOUT}.bam ${FILES_LIST} -@ ${THREADS}

# index
samtools index ${REPLICATE}.${LAYOUT}.bam -@ ${THREADS}

## bam to bedGraph, bedGraph to bigWig
#genomeCoverageBed -ibam ${REPLICATE}.${LAYOUT}.bam -bg -split > ${REPLICATE}.${LAYOUT}.bedGraph
#wigToBigWig ${REPLICATE}.${LAYOUT}.bedGraph ${CHR_LENGTH} ${REPLICATE}.${LAYOUT}.bw
#[ -f "${REPLICATE}.${LAYOUT}.bw" ] && rm ${REPLICATE}.${LAYOUT}.bedGraph
