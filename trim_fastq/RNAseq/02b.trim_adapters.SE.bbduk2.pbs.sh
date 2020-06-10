#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=30g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.trim_adapters.SE
#PBS -j oe
#PBS -J 0-8
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6
MEMORY=30g

IN_DIR=../Links
IN_SEQ=($(find $IN_DIR -maxdepth 1 -name "*.txt.gz" | grep -v "all"))
UNIQ_SEQ=($(printf "%s\n" "${IN_SEQ[@]%.txt.gz}" | sort -u))
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}

ADAPTER_DIR=/common/DB/genome_reference/other/adapters
ADAPTER="${ADAPTER_DIR}/Illumina_RNAseq.fa,${ADAPTER_DIR}/bbmap_adapters.fa"

BBDUK_PAR="overwrite=t \
ktrim=r \
k=23 \
rcomp=t \
mink=12 \
hdist=1 \
minoverlap=8 \
minlength=25 \
threads=$THREADS \
-Xmx$MEMORY"

# ----------------Commands------------------- #
# right trim
bbduk.sh in=${FILE}.txt.gz out=${BASE}.trim.txt.gz ref=$ADAPTER stats=${BASE}.stats $BBDUK_PAR 2> ${BASE}.trim.log
