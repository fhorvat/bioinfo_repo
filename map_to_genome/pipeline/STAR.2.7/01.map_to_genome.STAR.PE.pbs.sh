#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.map_to_genome.STAR
#PBS -l select=ncpus=12:mem=60g
#PBS -J 0-%N_SAMPLES
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
STAR_INDEX=%GENOME_DIR/STAR_index.2.7/%SJDB_OVERHANG
CHR_LENGTH=${STAR_INDEX}/chrNameLength.txt

INPUT_DIR=../../Raw/Links
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 \( -name "*.txt.gz" -not -name "*all*" \)))
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
STAR.2.7 --readFilesIn ${FILE}_1.txt.gz ${FILE}_2.txt.gz $STAR_PAR

# rename bam
mv ${BASE}.Aligned.sortedByCoord.out.bam ${BASE}.bam

# bam index
samtools index -@ $THREADS ${BASE}.bam
