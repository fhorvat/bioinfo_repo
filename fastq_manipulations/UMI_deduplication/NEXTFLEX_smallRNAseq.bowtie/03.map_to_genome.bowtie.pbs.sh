#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=18:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03.map_to_genome.bowtie
##PBS -J 0-7
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12

IN_DIR=.
IN_SEQ=($(find $IN_DIR -maxdepth 1 \( -name "*.no_space.fastq" -and -name "s_5oocytes_Mov10l1_KO_So869_5xoo_r1.SE.no_space.fastq" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.no_space.fastq}

IN_GENOME=/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1
IN_INDEX=${IN_GENOME}/bowtie_index/MesAur1.0

# ----------------Commands------------------- #
# bowtie allignment
bowtie --fullref -l 15 -n 2 -S -p 6 ${IN_INDEX} ${FILE} | samtools view -@ 6 -Sb - | samtools sort -@ 6 -o ${BASE}.sorted.bam -

# index 
samtools index ${BASE}.sorted.bam

