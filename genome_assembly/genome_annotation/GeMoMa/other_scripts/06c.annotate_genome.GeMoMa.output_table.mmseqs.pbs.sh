#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=5g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.06c.annotate_genome.GeMoMa
#PBS -J 0-1
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1
MEMORY=5G

INPUT_DIR=./mmseqs_database

IN_PROTEINS=($(find $INPUT_DIR -name "GCF_000*_genomic"))
FILE_PROTEINS=${IN_PROTEINS[$PBS_ARRAY_INDEX]}
BASE_PROTEINS=${FILE_PROTEINS#${INPUT_DIR}/}

IN_GENOME=($(find $INPUT_DIR -name "hamster.sequel.draft-20200302.arrow"))
FILE_GENOME=${IN_GENOME[0]}
BASE_GENOME=${FILE_GENOME#${INPUT_DIR}/}

# ----------------Commands------------------- #
# convert alignment to human-readable format 
mmseqs convertalis ${FILE_PROTEINS} ${FILE_GENOME} ./mmseqs.${BASE_PROTEINS}/mmseqs.${BASE_PROTEINS} ./mmseqs.${BASE_PROTEINS}/mmseqs.${BASE_PROTEINS}.tab \
--format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,target,raw,nident,empty,"
