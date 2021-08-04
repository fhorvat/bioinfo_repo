#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=40g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.sample_fastq
#PBS -j oe
#PBS -J 0-2
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# threads 
THREADS=6
MEMORY=40g

# input 
IN_DIR=..
IN_SEQ=($(find ${IN_DIR} -maxdepth 1 \( -name "*.PE_1.txt.gz" -not -name "*sample*" \)))
UNIQ_SEQ=($(printf "%s\n" "${IN_SEQ[@]%.PE_*.txt.gz}" | sort -u))
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.txt.gz}

# how many reads?
N=1000000

# ----------------Commands------------------- #
# get random N files
reformat.sh -Xmx$MEMORY in=${FILE}.PE_1.txt.gz in2=${FILE}.PE_2.txt.gz \
out=${BASE}.sample.PE_1.txt.gz out2=${BASE}.sample.PE_2.txt.gz \
samplereadstarget=$N sampleseed=13
