#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.map_to_genome.STAR
#PBS -l select=ncpus=10:mem=60g
#PBS -J 0-17
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=10
STAR_INDEX=/common/DB/genome_reference/cow/bosTau9.ARS-UCD1.2.GCA_002263795.2/STAR_index/sjdbOverhang_249
CHR_LENGTH=${STAR_INDEX}/chrNameLength.txt

INPUT_DIR=../../Raw/Links
IN_SEQ=(`find $INPUT_DIR -name "*txt.gz"`)
UNIQ_SEQ=(`printf "%s\n" "${IN_SEQ[@]%_[1-2].txt.gz}" | sort -u`)
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.txt.gz}

STAR_PAR="--genomeDir $STAR_INDEX \
--runThreadN $THREADS \
--genomeLoad NoSharedMemory \
--limitBAMsortRAM  20000000000 \
--readFilesCommand unpigz -c \
--outFileNamePrefix ${BASE}. \
--outSAMtype BAM SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFilterMultimapNmax 99999 \
--outFilterMultimapScoreRange 0 \
--outFilterMismatchNoverLmax 0.05 \
--sjdbScore 2"

# ----------------Commands------------------- #
# mapping
STAR --readFilesIn ${FILE} $STAR_PAR

# rename bam
mv ${BASE}.Aligned.sortedByCoord.out.bam ${BASE}.bam

# bam index
sambamba index ${BASE}.bam

# bam to bedGraph, bedGraph to bigWig
genomeCoverageBed -ibam ${BASE}.bam -bg -split -g $CHR_LENGTH > ${BASE}.bedGraph
wigToBigWig ${BASE}.bedGraph $CHR_LENGTH ${BASE}.bw
[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph
