#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=40g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.correct_kmers.PE
#PBS -j oe
#PBS -J 0-1
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6
MEMORY=40g

INPUT_DIR=../../Raw/Links
IN_SEQ=(`ls ${INPUT_DIR}/*.txt.gz`)
UNIQ_SEQ=(`printf "%s\n" "${IN_SEQ[@]%_*txt.gz}" | sort -u`)
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}

# ----------------Commands------------------- #
# unzip the files
#unpigz -p $THREADS -c ${FILE}_1.txt.gz > ${BASE}_1.fastq 
#unpigz -p $THREADS -c ${FILE}_2.txt.gz > ${BASE}_2.fastq

# correct erroneous k-mers
#perl ~/bin/run_rcorrector.pl -1 ${BASE}_1.fastq -2 ${BASE}_2.fastq -t $THREADS

# filter uncorrectable read pairs
python /common/WORK/fhorvat/programi/Rcorrector/FilterUncorrectabledPEfastq.py -1 ${BASE}_1.cor.fq -2 ${BASE}_2.cor.fq -o fixed

# remove old fastq
[ -f "fixed_${BASE}_1.cor.fq" ] && rm ${BASE}_1.cor.fq
[ -f "fixed_${BASE}_2.cor.fq" ] && rm ${BASE}_2.cor.fq

# rename
mv fixed_${BASE}_1.cor.fq ${BASE}_1.cor.fq
mv fixed_${BASE}_2.cor.fq ${BASE}_2.cor.fq

# zip
pigz -p $THREADS ${BASE}_1.cor.fq
pigz -p $THREADS ${BASE}_2.cor.fq

# rename
mv ${BASE}_1.cor.fq.gz ${BASE}_1.txt.gz
mv ${BASE}_2.cor.fq.gz ${BASE}_2.txt.gz

# remove original unzipped fastq
[ -f "${BASE}_1.txt.gz" ] && rm ${BASE}_1.fastq
[ -f "${BASE}_2.txt.gz" ] && rm ${BASE}_2.fastq
