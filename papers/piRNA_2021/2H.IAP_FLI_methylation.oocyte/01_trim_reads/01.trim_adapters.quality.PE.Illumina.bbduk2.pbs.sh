#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=40g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.trim_adapters.PE
#PBS -j oe
#PBS -J 0-2
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=60g

IN_DIR=../Links
IN_SEQ=($(find $IN_DIR -maxdepth 1 -name "*.txt.gz" | grep -v "all"))
UNIQ_SEQ=($(printf "%s\n" "${IN_SEQ[@]%_*txt.gz}" | sort -u))
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}

BBDUK_PAR="literal=AGATCGGAAGAGC \
overwrite=t \
ktrim=r \
k=12 \
rcomp=t \
mink=8 \
hdist=1 \
minoverlap=8 \
minlength=25 \
minlength=50 \
tbo \
threads=$THREADS \
-Xmx$MEMORY"

# ----------------Commands------------------- #
# right trim
bbduk.sh \
in1=${FILE}_1.txt.gz \
in2=${FILE}_2.txt.gz \
out1=${BASE}_1.trim.txt.gz \
out2=${BASE}_2.trim.txt.gz \
outs=${BASE}_s.trim.txt.gz \
stats=${BASE}.stats \
$BBDUK_PAR 2> ${BASE}.trim.log
