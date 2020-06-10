#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=50g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03.map_reads_to_transcriptome.minimap2
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=50G

BASE_DIR=`basename ${PWD/trinity*/}`

IN_DIR_READS=../../../Raw/${BASE_DIR}/Cleaned/trimmed
IN_SEQ_READS=(`ls ${IN_DIR_READS}/s_*all.PE_1.txt.gz`)
UNIQ_SEQ_READS=(`printf "%s\n" "${IN_SEQ_READS[@]%_*.txt.gz}" | sort -u`)
FILE_READS=${UNIQ_SEQ_READS[0]}

IN_DIR_TRANS=.
IN_SEQ_TRANS=(`ls ${IN_DIR_TRANS}/*.fasta`)
FILE_TRANS=${IN_SEQ_TRANS[0]}

BASE=${FILE_TRANS#${IN_DIR_TRANS}/}
BASE=${BASE/.fasta/.r2t}

MINIMAP2_PAR="-x sr \
--secondary=no \
-a \
-N 0 \
-t $THREADS"

# ----------------Commands------------------- #
# align reads to transcriptome
minimap2 $MINIMAP2_PAR ${FILE_TRANS} ${FILE_READS}_1.txt.gz ${FILE_READS}_2.txt.gz | samtools view -F 4 -F 2048 -@$THREADS -Sb -o ${BASE}.bam

# count mapped reads
samtools view ${BASE}.bam -@$THREADS | awk '{print $1}' | sort | uniq | wc -l > ${BASE}.read_counts.txt

# count all reads
source /common/WORK/fhorvat/Projekti/Svoboda/scripts/useful_functions.sh
countFastq ${FILE_READS}_1.txt.gz >> ${BASE}.read_counts.txt
