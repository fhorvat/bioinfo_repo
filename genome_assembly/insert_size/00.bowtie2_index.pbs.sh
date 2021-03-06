#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=4:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.00.bowtie2_index
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script parameters
THREADS=4

# set input variables
IN_GENOME=(`find . -name "*.fasta"`)
OUT_PREFIX=bowtie_index

# ----------------Commands------------------- #
# generate index
bowtie2-build --threads $THREADS $IN_GENOME $OUT_PREFIX

# create link to fasta 
ln -s $IN_GENOME ./${OUT_PREFIX}.fasta
