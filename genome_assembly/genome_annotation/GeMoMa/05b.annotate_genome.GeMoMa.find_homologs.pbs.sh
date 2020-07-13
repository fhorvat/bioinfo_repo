#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=18:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.05b.annotate_genome.GeMoMa
#PBS -J 0-1
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=18
MEMORY=20G

INPUT_DIR=.
IN_FASTA=($(find $INPUT_DIR -name "cds-parts.fasta"))
FILE_FASTA=${IN_FASTA[$PBS_ARRAY_INDEX]}
BASE_FASTA=${FILE_FASTA#${INPUT_DIR}/}
BASE_FASTA=${BASE_FASTA%/*}
BASE_FASTA=${BASE_FASTA#CDS.}

INPUT_DIR=../..
IN_GENOME=($(find $INPUT_DIR -name "*.fasta"))
FILE_GENOME=${IN_GENOME[0]}
BASE_GENOME=${FILE_GENOME#${INPUT_DIR}/}
BASE_GENOME=${BASE_GENOME%.fasta}

DATABASE_NAME=./tblastn_database.hamster.sequel.draft-20200302.arrow/${BASE_GENOME}

# ----------------Commands------------------- #
# create output dir
mkdir tblastn.${BASE_FASTA}

# extract introns and coverage from bam file
tblastn -outfmt "6 std sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles" \
-query ${FILE_FASTA} \
-db ${DATABASE_NAME} \
-out tblastn.${BASE_FASTA}/tblastn_results.txt \
-num_threads ${THREADS}

# sort by first column
sort -k1,1 tblastn.${BASE_FASTA}/tblastn_results.txt > tblastn.${BASE_FASTA}/tblastn_results.sorted.txt
