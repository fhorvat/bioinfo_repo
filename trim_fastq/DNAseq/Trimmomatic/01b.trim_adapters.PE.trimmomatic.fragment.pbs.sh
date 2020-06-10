#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=30g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01b.trim_adapters.PE
#PBS -j oe
#PBS -J 0-2
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6
MEMORY=30g

IN_DIR=../Links
IN_SEQ=($(find $IN_DIR -maxdepth 1 -name "s_mesAur_fragment*.txt.gz" | grep -v "all"))
UNIQ_SEQ=($(printf "%s\n" "${IN_SEQ[@]%_*txt.gz}" | sort -u))
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}

ADAPTERS=($(find . -name "Illumina_adapters.fasta"))
ADAPTER=${ADAPTERS[0]}

SCRIPT=/common/WORK/fhorvat/programi/Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar

# ----------------Commands------------------- #
# trim
java -jar $SCRIPT \
PE \
-threads ${THREADS} \
${FILE}_1.txt.gz \
${FILE}_2.txt.gz \
${BASE}_1.trim.txt.gz \
${BASE}_1.trim.un.txt.gz \
${BASE}_2.trim.txt.gz \
${BASE}_2.trim.un.txt.gz \
ILLUMINACLIP:${ADAPTER}:2:30:7:5:true \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:40 2> ${BASE}.log
