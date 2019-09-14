#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02a.merge_unmapped.map
#PBS -l select=ncpus=10:mem=60g
#PBS -J 0-%N_SAMPLES
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------MERGE------------------- #
INPUT_DIR=.
IN_SEQ=(`find $INPUT_DIR -name "*genome.Unmapped.out.mate1" | sed -e 's/.mate1//g'`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.Unmapped.out}.merged

# merge
bbmerge.sh in1=${FILE}.mate1 in2=${FILE}.mate2 out=${BASE}.fastq

# ----------------MAP MERGED------------------- #
# variables
THREADS=10
STAR_INDEX=%GENOME_DIR/STAR_index/%SJDB_OVERHANG
CHR_LENGTH=${STAR_INDEX}/chrNameLength.txt

FILE=${BASE}.fastq
BASE=${FILE%.fastq}

STAR_PAR="--genomeDir $STAR_INDEX \
--runThreadN $THREADS \
--genomeLoad NoSharedMemory \
--limitBAMsortRAM  20000000000 \
--outFileNamePrefix ${BASE}. \
--outSAMtype BAM SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFilterMultimapNmax 99999 \
--outFilterMultimapScoreRange 0 \
--outFilterMismatchNoverLmax 0.2 \
--sjdbScore 2"

# mapping
STAR --readFilesIn ${FILE} $STAR_PAR

# rename
mv ${BASE}.Aligned.sortedByCoord.out.bam ${BASE}.bam

