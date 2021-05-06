#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=60g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.sample_fastq
#PBS -j oe
#PBS -J 0-2
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# threads and memory
THREADS=12
HALF_THREADS=$(($THREADS / 2))

# input 
IN_DIR=.
IN_SEQ=($(find ${IN_DIR} \( -name "*.PE_1.txt.gz" -not -name "*sample*" \)))
UNIQ_SEQ=($(printf "%s\n" "${IN_SEQ[@]%.PE_*.txt.gz}" | sort -u))
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.txt.gz}

# how many reads?
N=100

# ----------------Commands------------------- #
# get random N files
# by Pierre Lindenbaum (https://www.biostars.org/p/6544/#6562)
paste <(unpigz -p ${HALF_THREADS} -c ${FILE}.PE_1.txt.gz) <(unpigz -p ${HALF_THREADS} -c ${FILE}.PE_2.txt.gz) |\
awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t\t");} }' |\
shuf  |\
head -n $N |\
sed 's/\t\t/\n/g' |\
awk -v NAME=$BASE -F '\t' '{print $1 > NAME".sample.PE_1.txt"; print $2 > NAME".sample.PE_2.txt"}'

## how it works, line-by-line:
# merge the two fastqs
# merge by group of 4 lines
# shuffle
# only 10 records
# restore the delimiters
# split in two files

# zip sampled fastq
pigz -p $THREADS ${BASE}.sample.PE_*.txt
