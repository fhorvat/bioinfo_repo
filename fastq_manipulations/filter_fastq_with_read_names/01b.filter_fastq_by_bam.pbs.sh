#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.00b.filter_fastq_by_read_name
#PBS -j oe
##PBS -J 0-3
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1
MEMORY=20G

# input bam file
IN_DIR=/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.smallRNAseq/Data/Mapped/STAR_Siomi.multimappers/6_filter_18to32nt
IN_DIR=.
IN_SEQ=($(find $IN_DIR -maxdepth 1 \( -name "*.bam" \)))
FILE_BAM=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE_BAM#${IN_DIR}/}
BASE=${BASE%.18to32nt.bam}

# input fastq file
IN_DIR=..
IN_SEQ=($(find ${IN_DIR} \( -name "*.atrim.txt.gz" -and -name "*${BASE}*" \)))
FILE_FASTQ=${IN_SEQ[0]}

# input read list
READ_LIST=$(find $IN_DIR -name "${BASE}.read_list.txt")

# ----------------Commands------------------- #
# filter fastq by read names
filterbyname.sh -Xmx$MEMORY in=${FILE_FASTQ} out=${BASE}.fastq names=<(samtools view -F 256 ${FILE_BAM} | awk '{print $1}') include=t
