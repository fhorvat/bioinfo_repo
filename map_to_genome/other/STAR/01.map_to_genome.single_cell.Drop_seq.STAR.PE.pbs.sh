#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.map_to_genome.STAR
#PBS -l select=ncpus=12:mem=60g
##PBS -J 0-48
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
STAR_INDEX=/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/STAR_index.2.7/sjdbOverhang_100
CHR_LENGTH=${STAR_INDEX}/chrNameLength.txt

INPUT_DIR=../../Raw/Links
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 \( -name "*.txt.gz" -not -name "*all*" \)))
UNIQ_SEQ=(`printf "%s\n" "${IN_SEQ[@]%_[1-2].txt.gz}" | sort -u`)
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}

#IN_READS_1=$(find $INPUT_DIR -maxdepth 1 \( -name "s_INT1*.PE_1.txt.gz" -not -name "*all*" \) | paste -sd"," -)
#IN_READS_2=$(find $INPUT_DIR -maxdepth 1 \( -name "s_INT1*.PE_2.txt.gz" -not -name "*all*" \) | paste -sd"," -)

STAR_PAR="--genomeDir $STAR_INDEX \
--runThreadN $THREADS \
--genomeLoad NoSharedMemory \
--limitBAMsortRAM  20000000000 \
--readFilesCommand unpigz -c \
--outSAMtype BAM SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFilterMultimapNmax 10 \
--outFilterMultimapScoreRange 0 \
--outFilterMismatchNoverLmax 0.05 \
--sjdbScore 2 \
--soloType CB_UMI_Simple \
--soloCBwhitelist None \
--soloCBstart 1 \
--soloCBlen 12 \
--soloUMIstart 13 \
--soloUMIlen 8 \
--soloBarcodeReadLength 26 \
--soloStrand Unstranded \
--soloFeatures Gene GeneFull SJ Velocyto \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM"

# ----------------Commands------------------- #
# mapping
STAR.2.7 --readFilesIn ${FILE}_1.txt.gz ${FILE}_2.txt.gz $STAR_PAR

# rename bam
#mv ${BASE}.Aligned.sortedByCoord.out.bam ${BASE}.bam

# bam index
#samtools index -@ $THREADS ${BASE}.bam
