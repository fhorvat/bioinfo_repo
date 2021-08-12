#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=60g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.annotate_genome.Liftoff
#PBS -J 0-4
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=60G

# target genome
INPUT_DIR=../..
IN_GENOME=($(find $INPUT_DIR -maxdepth 1 -name "*.fa"))
FILE_GENOME=${IN_GENOME[0]}
BASE_GENOME=${FILE_GENOME#${INPUT_DIR}/}
BASE_GENOME=${BASE_GENOME%.fa}

# query genome
INPUT_DIR=../gastropod_annotation
IN_QUERY=($(find $INPUT_DIR -maxdepth 1 -name "*genomic.fna"))
FILE_QUERY=${IN_QUERY[$PBS_ARRAY_INDEX]}
BASE_QUERY=${FILE_QUERY#${INPUT_DIR}/}
BASE_QUERY=${BASE_QUERY%_genomic.fna}

# query .gtf
FILE_GTF=${FILE_QUERY/.fna/.gff}

# script
SCRIPT=/common/WORK/fhorvat/programi/python/packages/bin/liftoff

# ----------------Commands------------------- #
# map annotation to another genome
$SCRIPT \
-g ${FILE_GTF} \
-p ${THREADS} \
-o ./${BASE_GENOME}.${BASE_QUERY}.liftoff.gff \
-dir ./${BASE_GENOME}.${BASE_QUERY}.tmp \
-u ./${BASE_GENOME}.${BASE_QUERY}.unmapped_features.txt \
${FILE_GENOME} \
${FILE_QUERY}
