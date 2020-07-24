#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=15g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.map_reads_to_sequence.bowtie
#PBS -J 0-8
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6
MEMORY=15G

IN_DIR=../../Raw/Cleaned
IN_SEQ=(${IN_DIR}/*.txt.gz)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.txt.gz}

INDEX=/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/bowtie2_index/mm10

BOWTIE2_PAR="--local \
-p $THREADS \
-q \
-x $INDEX"

# ----------------Commands------------------- #
# map with bowtie
bowtie2 $BOWTIE2_PAR \
-U ${FILE} 2> \
${BASE}.stats.txt | \
samtools view -@$THREADS -Sb - |
samtools sort - -T ${BASE} -o ${BASE}.bam
