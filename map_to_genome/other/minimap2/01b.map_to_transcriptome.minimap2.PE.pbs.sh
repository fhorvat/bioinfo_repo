#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=50g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.align.minimap2.PE
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=50G

BASE_DIR=`basename ${PWD/trinity*/}`
IN_DIR=../../../../Raw/$BASE_DIR/Cleaned/extended_merge.bbnorm_normalized
IN_SEQ=(`ls ${IN_DIR}/s_*all.PE_1.txt.gz`)
UNIQ_SEQ=(`printf "%s\n" "${IN_SEQ[@]%_*.txt.gz}" | sort -u`)
FILE=${UNIQ_SEQ[0]}
BASE=${FILE#${IN_DIR}/}

TRANSCRIPTOME_FASTA=`ls ../*fasta`

MINIMAP2_PAR="-x sr \
--secondary=no \
-a \
-N 0 \
-t $THREADS"

# ----------------Commands------------------- #
# align reads to transcriptome
minimap2 $MINIMAP2_PAR $TRANSCRIPTOME_FASTA ${FILE}_1.txt.gz ${FILE}_2.txt.gz | samtools view -F 4 -F 2048 -@$THREADS -Sb -o ${BASE}.bam

# count mapped reads
samtools view ${BASE}.bam -@$THREADS | awk '{print $1}' | sort | uniq | wc -l > ${BASE}.read_counts.txt

# count all reads
countFastq ${FILE} >> ${BASE}.read_counts.txt
