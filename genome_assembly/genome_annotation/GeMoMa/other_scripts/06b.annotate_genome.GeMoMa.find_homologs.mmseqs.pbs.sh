#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=18:mem=60g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.06b.annotate_genome.GeMoMa
#PBS -J 0-1
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=18
MEMORY=20G

INPUT_DIR=./mmseqs_database

IN_PROTEINS=($(find $INPUT_DIR -name "GCF_000*_genomic"))
FILE_PROTEINS=${IN_PROTEINS[$PBS_ARRAY_INDEX]}
BASE_PROTEINS=${FILE_PROTEINS#${INPUT_DIR}/}

IN_GENOME=($(find $INPUT_DIR -name "hamster.sequel.draft-20200302.arrow"))
FILE_GENOME=${IN_GENOME[0]}
BASE_GENOME=${FILE_GENOME#${INPUT_DIR}/}

# ----------------Commands------------------- #
# create output dir
mkdir ./mmseqs.${BASE_PROTEINS}

# extract introns and coverage from bam file
mmseqs search ${FILE_PROTEINS} ${FILE_GENOME} ./mmseqs.${BASE_PROTEINS}/mmseqs.${BASE_PROTEINS} ./tmp --threads ${THREADS}

## sort by first column
#sort -k1,1 tblastn.${BASE_FASTA}/tblastn_results.txt > tblastn.${BASE_FASTA}/tblastn_results.sorted.txt
