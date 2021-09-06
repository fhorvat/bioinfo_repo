#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.map_to_genome.STAR
#PBS -l select=ncpus=24:mem=60g
#PBS -j oe
#PBS -J 0-3
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# override number of threads 
THREADS=24

# input
INPUT_DIR=..
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 \( -name "*.txt.gz" -not -name "*atrim.txt.gz" \)))
UNIQ_SEQ=($(printf "%s\n" "${IN_SEQ[@]%_[1-2].txt.gz}" | sort -u))
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.txt.gz}

STAR_INDEX=/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed/STAR_index.2.7/sjdbOverhang_100

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
STAR.2.7 --readFilesIn ${FILE} ${STAR_PAR}

# rename bam
mv ${BASE}.Aligned.sortedByCoord.out.bam ${BASE}.mapped.bam
