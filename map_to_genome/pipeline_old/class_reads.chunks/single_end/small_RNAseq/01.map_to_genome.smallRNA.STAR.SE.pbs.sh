#/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.map_to_genome.STAR
#PBS -l select=ncpus=10:mem=60g
#PBS -J 0-%N_SAMPLES
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=8
STAR_INDEX=%GENOME_DIR/STAR_index/%SJDB_OVERHANG
CHR_LENGTH=$STAR_INDEX/chrNameLength.txt

INPUT_DIR=../../Raw/Cleaned
IN_SEQ=(`ls ${INPUT_DIR}/*.txt.gz`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.txt.gz}.genome

STAR_PAR="--genomeDir $STAR_INDEX \
--runThreadN $THREADS \
--genomeLoad LoadAndRemove \
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
--outFilterMultimapNmax 99999 \
--outFilterMultimapScoreRange 0 \
--alignIntronMax 1 \
--alignSJDBoverhangMin 999999999999"

# ----------------Commands------------------- #
# mapping
STAR --readFilesIn ${FILE} $STAR_PAR

# bam index
samtools index ${BASE}.Aligned.sortedByCoord.out.bam

# bam to bedGraph, bedGraph to bigWig
genomeCoverageBed -ibam ${BASE}.Aligned.sortedByCoord.out.bam -bg -split -g $CHR_LENGTH > ${BASE}.bedGraph
wigToBigWig ${BASE}.bedGraph $CHR_LENGTH ${BASE}.bw
[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph
