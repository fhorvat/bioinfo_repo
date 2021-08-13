#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.00a.bam_read_list
#PBS -J 0-3
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
IN_DIR=/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.smallRNAseq/Data/Mapped/STAR_Siomi.multimappers/6_filter_18to32nt
IN_SEQ=($(find $IN_DIR -maxdepth 1 \( -name "*.bam" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.bam}

# ----------------Commands------------------- #
# get all unique reads from bam as a list
#sambamba view ${FILE} | awk '{ a[$1]++ } END { for (b in a) { print b } }' > ${BASE}.read_list.txt

# get all primary alignments read names from bam as a list
samtools view -F 256 ${FILE} | awk '{print $1}' > ${BASE}.read_list.txt
