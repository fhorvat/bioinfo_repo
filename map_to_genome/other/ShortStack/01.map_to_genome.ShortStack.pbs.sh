#/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.map_to_genome.ShortStack
#PBS -l select=ncpus=10:mem=60g
#PBS -J 0-1
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=10
GENOME=./mm10.fa

INPUT_DIR=./Links
IN_SEQ=(`ls ${INPUT_DIR}/*.fastq.gz`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fastq.gz}

# ----------------Commands------------------- #
# mapping
ShortStack --outdir $BASE --bowtie_cores $THREADS --sort_mem 40G --nohp --pad 0 --mmap u --readfile $FILE --genomefile $GENOME

# bam index
#samtools index ${BASE}.Aligned.sortedByCoord.out.bam

# bam to bedGraph, bedGraph to bigWig
#genomeCoverageBed -ibam ${BASE}.Aligned.sortedByCoord.out.bam -bg -split -g $CHR_LENGTH > ${BASE}.bedGraph
#wigToBigWig ${BASE}.bedGraph $CHR_LENGTH ${BASE}.bw
#[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph

