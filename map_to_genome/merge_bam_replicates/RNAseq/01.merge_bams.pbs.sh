#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.merge_bams
#PBS -l select=ncpus=12:mem=30g
#PBS -J 0-%N_REPLICATES
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------MAP MERGED------------------- #
# source
source ./000.load_variables.sh

# override number of threads
THREADS=12

# input
REPLICATE=${REPLICATES[$PBS_ARRAY_INDEX]}
FILES=($(find ${INPUT_DIR} -maxdepth 1 \( -name "${REPLICATE}*r*.*E.*bam" \)))
FILES=$(printf "%s " ${FILES[@]})

# just in case
echo $FILES

# ----------------Commands------------------- #
# merge
samtools merge ${REPLICATE}.bam ${FILES} -@ ${THREADS}

# index
samtools index ${REPLICATE}.bam -@ ${THREADS}

## bam to bedGraph, bedGraph to bigWig
#genomeCoverageBed -ibam ${STAGE}.bam -bg -split > ${REPLICATE}.bedGraph
#wigToBigWig ${REPLICATE}.bedGraph ${CHR_LENGTH} ${REPLICATE}.bw
#[ -f "${REPLICATE}.bw" ] && rm ${REPLICATE}.bedGraph
