#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.merge_unmapped.map
#PBS -l select=ncpus=10:mem=60g
#PBS -J 0-%N_SAMPLES
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------MERGE------------------- #
# variables
INPUT_DIR=.
IN_SEQ=(`ls $INPUT_DIR/*.mate* | grep -v "45S" | grep -v "merged"`)
UNIQ_SEQ=(`printf "%s\n" "${IN_SEQ[@]%.Unmapped.out.mate*}" | sort -u`)
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}.merged

# merging
echo "--- a) merge unmapped ---"
bbmerge.sh in1=${FILE}.Unmapped.out.mate1 in2=${FILE}.Unmapped.out.mate2 out=${BASE}.fastq

# ----------------MAP MERGED------------------- #
# variables
THREADS=8
STAR_INDEX=%GENOME_DIR/STAR_index/%SJDB_OVERHANG

FILE=${BASE}.fastq
BASE=${FILE%.fastq}
STAR_PAR="--genomeDir $STAR_INDEX \
--runThreadN $THREADS \
--genomeLoad LoadAndRemove \
--limitBAMsortRAM  20000000000 \
--outFileNamePrefix ${BASE}. \
--outSAMtype BAM SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFilterMultimapNmax 99999 \
--outFilterMultimapScoreRange 0 \
--outFilterMismatchNoverLmax 0.2 \
--sjdbScore 2"

# mapping
echo "--- b) map merged unmapped ---"
STAR --readFilesIn ${FILE} $STAR_PAR
