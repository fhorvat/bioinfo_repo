#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=10g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.bam_to_fastq
#PBS -j oe
#PBS -J 0-2
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1
MEMORY=10G

IN_DIR=.
IN_SEQ=($(find ${IN_DIR} -maxdepth 1 \( -name "*.bam" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.bam}
BASE=${BASE/_bismark_bt2_pe.deduplicated/}

# ----------------Commands------------------- #
# bam to fastq
samtools fastq \
-0 ${BASE}_0.txt \
-1 ${BASE}_1.txt \
-2 ${BASE}_2.txt \
-s ${BASE}_s.txt \
$FILE

# gzip
for file in `find . -name "${BASE}_*.txt"`; do pigz -p 1 $file; done
