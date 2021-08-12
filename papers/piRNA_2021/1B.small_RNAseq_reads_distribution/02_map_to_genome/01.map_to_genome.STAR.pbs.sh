#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.map_to_genome.STAR
#PBS -l select=ncpus=12:mem=60g
#PBS -j oe
#PBS -J 0-4
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
STAR_PAR="--genomeDir $STAR_INDEX \
--runThreadN $THREADS \
--genomeLoad NoSharedMemory \
--limitBAMsortRAM  20000000000 \
--readFilesCommand unpigz -c \
--outFileNamePrefix ${BASE}. \
--outSAMtype BAM SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFilterMismatchNmax 1 \
--outFilterMismatchNoverLmax 1 \
--outFilterMismatchNoverReadLmax 1 \
--outFilterMatchNmin 16 \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
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
#--alignEndsType EndToEnd \
--alignIntronMax 1 \
--alignSJDBoverhangMin 999999999999"

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
