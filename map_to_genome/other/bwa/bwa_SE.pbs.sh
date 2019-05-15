#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=80g
#PBS -J 0-9
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N bwa
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=4

REF_DIR=/common/WORK/fhorvat/reference/mouse/mm10
REF_GENOME=${REF_DIR}/genome_indexes/bwa/mm10.fa
CHR_SIZE=${REF_DIR}/chrNameLength.txt

INPUT_DIR=/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/Srinivasan_2016_brain/Data/Raw/Links
IN_SEQ=($INPUT_DIR/*_neuron*.fastq.gz)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE##*$INPUT_DIR/}
BASE=${BASE%.fastq.gz}

# ----------------Commands------------------- #
# mapping 
bwa mem -t $THREADS $REF_GENOME $FILE > ${BASE}.sam

# sam to bam
samtools view ${BASE}.sam -@ $THREADS -Sb > ${BASE}.bam
samtools sort -@ $THREADS ${BASE}.bam ${BASE}_sorted
samtools index ${BASE}_sorted.bam
[ -f "${BASE}_sorted.bam" ] && rm -f ${BASE}.bam
[ -f "${BASE}_sorted.bam" ] && rm -f ${BASE}.sam

# bam to bedGraph,  bedGraph to bigWig
genomeCoverageBed -ibam ${BASE}_sorted.bam -bg -split -g $CHR_SIZE > ${BASE}.bedGraph
wigToBigWig ${BASE}.bedGraph $REF_DIR/chrNameLength.txt ${BASE}.bw
[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph
