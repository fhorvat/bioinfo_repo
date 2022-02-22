#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.map_to_genome.STAR
#PBS -l select=ncpus=12:mem=60g
#PBS -J 0-7
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# source
source ./000.load_variables.sh

# override number of threads 
THREADS=12

# input
UNIQ_SEQ=($(printf "%s\n" "${IN_SEQ[@]%_[1-2].txt.gz}" | sort -u))
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.txt.gz}

# mapping parameters
if [ ${ALL_MULTIMAPPERS} == FALSE ]
then
        STAR_PAR="--genomeDir ${STAR_INDEX} \
        --runThreadN ${THREADS} \
        --genomeLoad NoSharedMemory \
        --limitBAMsortRAM  20000000000 \
        --readFilesCommand unpigz -c \
        --outFileNamePrefix ${BASE}. \
        --outSAMtype BAM SortedByCoordinate \
        --outReadsUnmapped Fastx \
        --outFilterMultimapNmax 20 \
        --outFilterMultimapScoreRange 0 \
        --outFilterMismatchNoverLmax 0.05 \
        --sjdbScore 2"
else
        STAR_PAR="--genomeDir ${STAR_INDEX} \
        --runThreadN ${THREADS} \
        --genomeLoad NoSharedMemory \
        --limitBAMsortRAM  20000000000 \
        --readFilesCommand unpigz -c \
        --outFileNamePrefix ${BASE}. \
        --outSAMtype BAM SortedByCoordinate \
        --outReadsUnmapped Fastx \
        --outFilterMultimapNmax 5000 \
        --winAnchorMultimapNmax 5000 \
        --seedSearchStartLmax 30 \
        --alignTranscriptsPerReadNmax 30000 \
        --alignWindowsPerReadNmax 30000 \
        --alignTranscriptsPerWindowNmax 300 \
        --seedPerReadNmax 3000 \
        --seedPerWindowNmax 300 \
        --seedNoneLociPerWindow 1000 \
        --outFilterMultimapScoreRange 0 \
        --outFilterMismatchNoverLmax 0.05 \
        #--alignEndsType EndToEnd \
        --sjdbScore 2"
fi

# ----------------Commands------------------- #
# mapping
if [ ${SINGLE_END} == TRUE ]
then
   STAR.2.7 --readFilesIn ${FILE} ${STAR_PAR}
else
   STAR.2.7 --readFilesIn ${FILE}_1.txt.gz ${FILE}_2.txt.gz ${STAR_PAR}
fi

# rename bam
mv ${BASE}.Aligned.sortedByCoord.out.bam ${BASE}.bam

# bam index
samtools index -@ ${THREADS} ${BASE}.bam
