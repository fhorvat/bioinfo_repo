#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.merge_bams
#PBS -l select=ncpus=6:mem=30g
#PBS -J 0-%N_JOBS
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------MAP MERGED------------------- #
# source
source ./000.load_variables.sh

# override number of threads
THREADS=12

# replicate index
REPLICATE_INDEX=$((${PBS_ARRAY_INDEX} / ${N_READ_LENGTHS}))
REPLICATE=${REPLICATES[${REPLICATE_INDEX}]}

# read length index
READ_LENGTH_INDEX=$((${PBS_ARRAY_INDEX} % ${N_READ_LENGTHS}))
READ_LENGTH=${READ_LENGTHS[${READ_LENGTH_INDEX}]}

# print 
echo "replicate index:" ${REPLICATE_INDEX}", replicate:" ${REPLICATE}
echo "read length index:" ${READ_LENGTH_INDEX}", read length:" ${READ_LENGTH}

# ----------------Commands------------------- #
# loop through different read lengths and merge
if [ ${READ_LENGTH} == "all" ]
then
	FILES_LIST=($(find ${INPUT_DIR} -maxdepth 1 \( -name "${REPLICATE}*.[S,P]E.bam" -not -name "*bai" \)))
        FILES_LIST=$(printf "%s " ${FILES_LIST[@]})
        FILE_NAME=${REPLICATE}.${LAYOUT}
else
	FILES_LIST=($(find ${INPUT_DIR} -maxdepth 1 \( -name "${REPLICATE}*.${READ_LENGTH}.bam" -not -name "*bai" \)))
        FILES_LIST=$(printf "%s " ${FILES_LIST[@]})
        FILE_NAME=${REPLICATE}.${LAYOUT}.${READ_LENGTH}
fi
        
# just in case
printf "%s\n" ${FILES_LIST}

# merge
samtools merge ${FILE_NAME}.bam ${FILES_LIST} -@ ${THREADS}

# index
samtools index ${FILE_NAME}.bam -@ ${THREADS}

