#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.05a.annotate_genome.GeMoMa
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1
MEMORY=20G

INPUT_DIR=../../
IN_GENOME=($(find $INPUT_DIR -name "*.fasta"))
FILE_GENOME=${IN_GENOME[0]}
BASE_GENOME=${FILE_GENOME#${INPUT_DIR}}
BASE_GENOME=${BASE_GENOME%.fasta}

# ----------------Commands------------------- #
# create output dir
mkdir tblastn_database.${BASE_GENOME}

# extract introns and coverage from bam file
makeblastdb \
-in ${FILE_GENOME} \
-dbtype nucl \
-title ${BASE_GENOME} \
-out ./tblastn_database.${BASE_GENOME}
