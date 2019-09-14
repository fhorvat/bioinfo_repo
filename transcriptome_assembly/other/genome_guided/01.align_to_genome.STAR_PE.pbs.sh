#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.align_to_genome.STAR
#PBS -l select=ncpus=10:mem=60g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=8
STAR_INDEX=/common/WORK/fhorvat/reference/mouse/mm10/STAR_index/sjdbOverhang_249
CHR_LENGTH=${STAR_INDEX}/chrNameLength.txt

INPUT_DIR=..
BASE="pwd"

STAR_PAR="--genomeDir $STAR_INDEX \
--runThreadN $THREADS \
--genomeLoad LoadAndRemove \
--limitBAMsortRAM  20000000000 \
--readFilesCommand unpigz -c \
--outFileNamePrefix ${BASE}. \
--outSAMtype BAM SortedByCoordinate \
--outReadsUnmapped None \
--outFilterMultimapNmax 99999 \
--outFilterMismatchNoverLmax 0.2 \
--sjdbScore 2"

# ----------------Commands------------------- #
# concatenate fastq files to one file
cat ${INPUT_DIR}/*txt.gz > ${BASE}.merged.txt.gz

# mapping
STAR --readFilesIn ${BASE}.merged.txt.gz $STAR_PAR
