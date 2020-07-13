#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=10g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03.bam_to_fastq
#PBS -j oe
#PBS -J 0-1
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1
MEMORY=10G

IN_DIR=.
IN_SEQ=($(find ${IN_DIR} -maxdepth 1 \( -name "*.bam" -and -name "*merged*" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE%.merged.bam}

# ----------------Commands------------------- #
# bam to fastq
samtools fastq \
-0 ${BASE}_0.fastq \
-1 ${BASE}_1.fastq \
-2 ${BASE}_2.fastq \
-s ${BASE}_s.fastq \
$FILE
