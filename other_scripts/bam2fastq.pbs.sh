#!/bin/bash
# INFO: takes primary alignments from bam, outputs fastq and merges with unmapped from STAR
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.bam2fastq
#PBS -l select=ncpus=1:mem=20g
#PBS -J 0-15
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
INPUT_DIR=.
IN_SEQ=($INPUT_DIR/*.Aligned.sortedByCoord.out.bam)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.Aligned.sortedByCoord.out.bam}

# ----------------Commands------------------- #
# bam to fastq
samtools view -F 256 -b ${FILE} | bamToFastq -i - -fq ${BASE}.fastq

# merge with unmapped
cat ${BASE}.Unmapped.out.mate1 >> ${BASE}.fastq

# counts reads in fastq
echo ${BASE}.fastq `awk '{s++}END{print s/4}' ${BASE}.fastq` >> stats.txt

# gzip fastq
mv ${BASE}.fastq ${BASE}.txt
gzip ${BASE}.txt
