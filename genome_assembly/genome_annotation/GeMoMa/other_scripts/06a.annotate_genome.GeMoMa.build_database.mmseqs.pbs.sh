#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.06a.annotate_genome.GeMoMa
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1
MEMORY=20G

INPUT_DIR=../..
IN_GENOME=($(find $INPUT_DIR -name "*.fasta"))
FILE_GENOME=${IN_GENOME[0]}
BASE_GENOME=${FILE_GENOME#${INPUT_DIR}/}
BASE_GENOME=${BASE_GENOME%.fasta}

INPUT_DIR=.
IN_FASTA=($(find $INPUT_DIR -name "cds-parts.fasta"))

FILE_FASTA_1=${IN_FASTA[0]}
BASE_FASTA_1=${FILE_FASTA_1#${INPUT_DIR}/}
BASE_FASTA_1=${BASE_FASTA_1%/*}
BASE_FASTA_1=${BASE_FASTA_1#CDS.}

FILE_FASTA_2=${IN_FASTA[1]}
BASE_FASTA_2=${FILE_FASTA_2#${INPUT_DIR}/}
BASE_FASTA_2=${BASE_FASTA_2%/*}
BASE_FASTA_2=${BASE_FASTA_2#CDS.}

# ----------------Commands------------------- #
# create output dir
mkdir ./mmseqs_database

# create genome database
mmseqs createdb ${FILE_GENOME} ./mmseqs_database/${BASE_GENOME}

# create CDS databases
mmseqs createdb ${FILE_FASTA_1} ./mmseqs_database/${BASE_FASTA_1}
mmseqs createdb ${FILE_FASTA_2} ./mmseqs_database/${BASE_FASTA_2}
