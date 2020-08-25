#/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.map_to_genome.STAR
#PBS -l select=ncpus=12:mem=60g
#PBS -J 0-21
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
STAR_INDEX=/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/STAR_index/sjdbOverhang_249
CHR_LENGTH=$STAR_INDEX/chrNameLength.txt

INPUT_DIR=../../Raw/Cleaned
IN_SEQ=(${INPUT_DIR}/*.txt.gz)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.txt.gz}

STAR_PAR="--genomeDir $STAR_INDEX \
--runThreadN $THREADS \
#--genomeLoad LoadAndRemove \
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
--outFilterMultimapNmax 99999 \
--outFilterMultimapScoreRange 0 \
--alignIntronMax 1 \
--alignSJDBoverhangMin 999999999999"

# ----------------Commands------------------- #
# mapping
STAR --readFilesIn ${FILE} $STAR_PAR

# rename bam
mv ${BASE}.Aligned.sortedByCoord.out.bam ${BASE}.bam

# index
samtools index ${BASE}.bam
