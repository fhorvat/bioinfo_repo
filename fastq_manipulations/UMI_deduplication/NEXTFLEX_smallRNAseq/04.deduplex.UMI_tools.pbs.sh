#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.04.deduplex.umi_tools
##PBS -J 0-7
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12

IN_DIR=.
IN_SEQ=($(find $IN_DIR -maxdepth 1 \( -name "*.sorted.bam" -and -name "s_5oocytes_Mov10l1_KO_So869_5xoo_r1.SE.sorted.bam" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.sorted.bam}

# ----------------Commands------------------- #
# UMItools
umi_tools dedup --method directional -I ${BASE}.sorted.bam -S ${BASE}.deduped.bam

# convert deduped bam files to fastq files
bam2fastx -q -Q -A -o ${BASE}.deduped.fastq ${BASE}.deduped.bam
