#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=50g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.align.minimap2.SE
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=50G

BASE_DIR=`basename ${PWD/trinity*/}`
INPUT_DIR=../../../../Raw/$BASE_DIR/Cleaned/extended_merge.bbnorm_normalized
IN_SEQ=(`ls ${INPUT_DIR}/s_*all.PE.txt.gz`)
FILE=${IN_SEQ[0]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.PE.txt.gz}

TRANSCRIPTOME_FASTA=`ls ../*fasta`

MINIMAP2_PAR="-x sr \
--secondary=no \
-a \
-N 0 \
-t $THREADS"

# ----------------Commands------------------- #
# align reads to transcriptome
minimap2 $MINIMAP2_PAR $TRANSCRIPTOME_FASTA ${FILE} | samtools view -F 4 -F 2048 -@$THREADS -Sb -o ${BASE}.SE.bam

# split mapped and unmapped reads
#samtools view ${BASE}.bam -F 4 -F 2048 -@$THREADS -Sb -o ${BASE}.mapped.bam
#samtools view ${BASE}.bam -F 4 -@$THREADS -Sb -o ${BASE}.mapped.bam
#samtools view ${BASE}.bam -f 4 -@$THREADS -Sb -o ${BASE}.unmapped.bam

# count mapped reads
samtools view ${BASE}.SE.bam -@$THREADS | awk '{print $1}' | sort | uniq | wc -l > ${BASE}.SE.read_counts.txt

# count all reads
countFastq ${FILE} >> ${BASE}.SE.read_counts.txt
