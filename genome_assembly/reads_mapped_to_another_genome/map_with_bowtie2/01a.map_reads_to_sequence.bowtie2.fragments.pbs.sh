#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=15g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01a.map_reads_to_sequence.bowtie
#PBS -J 0-3
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=15G

IN_DIR=../../../../../../Raw/Cleaned
IN_SEQ=(${IN_DIR}/s_mesAur_fragment*.PE_1.txt.gz)
UNIQ_SEQ=($(printf "%s\n" "${IN_SEQ[@]%_*.txt.gz}" | sort -u))
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.txt.gz}

INSERT_LOW=50
INSERT_HIGH=300

INDEX=bowtie2_index

BOWTIE2_PAR="--fast-local \
-p $THREADS \
-q \
--fr \
-I $INSERT_LOW \
-X $INSERT_HIGH \
-x $INDEX \
-N 1 --ma 3 --mp 5 --rdg 4,2 --rfg 4,2 --score-min G,18,7"

# ----------------Commands------------------- #
# map with bowtie
bowtie2 $BOWTIE2_PAR \
-1 ${FILE}_1.txt.gz \
-2 ${FILE}_2.txt.gz \
-U ${FILE}_s.txt.gz 2> \
${BASE}.stats.txt | \
samtools view -@$THREADS -Sb -F 4 - |
samtools sort -@$THREADS - -T ${BASE} -o ${BASE}.bam 
