#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=40g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01a.trim_adapters.PE
#PBS -j oe
#PBS -J 0-8
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=60g

IN_DIR=../Links
IN_SEQ=($(find $IN_DIR -maxdepth 1 -name "*.txt.gz" | grep -v "all"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.txt.gz}

ADAPTERS=($(find . -name "Illumina_adapters.fasta"))
ADAPTER=${ADAPTERS[0]}

BBDUK_PAR="overwrite=t \
ktrim=r \
k=19 \
rcomp=t \
mink=8 \
hdist=1 \
minoverlap=8 \
minlength=25 \
qtrim=rl \
trimq=25 \
minlength=50 \
minavgquality=30 \
tbo \
threads=$THREADS \
-Xmx$MEMORY"

# ----------------Commands------------------- #
# right trim
bbduk.sh \
in=${FILE} \
out=${BASE}.txt.gz \
ref=$ADAPTER \
stats=${BASE}.stats \
$BBDUK_PAR 2> ${BASE}.log
