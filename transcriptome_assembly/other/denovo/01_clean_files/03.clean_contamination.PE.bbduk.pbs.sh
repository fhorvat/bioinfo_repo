#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=30g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03.clean_contamination.PE
#PBS -j oe
#PBS -J 0-2
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6
MEMORY=30g

REFERENCE=/common/DB/genome_reference/adapters/phix_adapters.fa

IN_DIR=.
IN_SEQ=(`ls ${IN_DIR}/*.trim.txt.gz`)
UNIQ_SEQ=(`printf "%s\n" "${IN_SEQ[@]%_*.trim.txt.gz}" | sort -u`)
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}

BBDUK_PAR="ref=$REFERENCE \
overwrite=t \
k=28 \
hdist=1 \
t=$THREADS \
-Xmx$MEMORY"

# ----------------Commands------------------- #
bbduk.sh \
in1=${FILE}_1.trim.txt.gz in2=${FILE}_2.trim.txt.gz \
out1=${BASE}_1.clean.trim.txt.gz out2=${BASE}_2.clean.trim.txt.gz \
outm1=${BASE}_1.dirty.trim.txt.gz outm2=${BASE}_2.dirty.trim.txt.gz \
stats=${BASE}.clean.trim.stats $BBDUK_PAR 2> ${BASE}.clean.trim.log
