#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03.library_sizes
#PBS -J 0-39
#PBS -l select=ncpus=1:mem=20g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
INPUT_DIR=.
IN_SEQ=($INPUT_DIR/*.bam)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}

# ----------------Commands------------------- #
# get histogram of CIGAR values of unique alignments
samtools view ${BASE}.bam | awk '!($1 in a){a[$1]; print $6}' | sort -n | uniq -c | sed -e 's/^ *//g;s/ /\'$'\t/g'  > library_hist.${BASE}.txt

# sum all alignments count to get total library size, add to one file
LIB_SIZE=`awk '{total+=$1}END{print total}' library_hist.${BASE}.txt`
echo -e $BASE'\t'$LIB_SIZE >> library_sizes.txt
