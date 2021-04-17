#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N STAR
#PBS -l select=ncpus=6:mem=40g
#PBS -J 0-4
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=4
REF_DIR=/common/WORK/fhorvat/reference/mouse/mm10
GENOME_DIR=${REF_DIR}/genome_indexes/STAR/sjdbOverhang_84
CHR_LENGTH=${GENOME_DIR}/chrNameLength.txt

INPUT_DIR=/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Data/Raw/Links
IN_SEQ=(`find $INPUT_DIR -name "WT*.fastq.gz" -o -name "Lnc5*.fastq.gz"`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE##*$INPUT_DIR/}
BASE=${BASE%.fastq.gz}

STAR_PAR="--genomeDir $GENOME_DIR \
--runThreadN $THREADS \
--genomeLoad LoadAndRemove \
--limitBAMsortRAM  20000000000 \
--readFilesCommand unpigz -c \
--outFileNamePrefix ${BASE}_"

# ----------------Commands------------------- #
# mapping
STAR --readFilesIn ${FILE} $STAR_PAR &> ${BASE}.log
