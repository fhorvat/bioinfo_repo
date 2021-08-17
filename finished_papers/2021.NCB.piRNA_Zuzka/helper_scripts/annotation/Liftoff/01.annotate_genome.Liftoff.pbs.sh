#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=60g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.annotate_genome.Liftoff
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=60G

# target genome
INPUT_DIR=/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed
IN_GENOME=($(find $INPUT_DIR -maxdepth 1 -name "*.fasta"))
FILE_GENOME=${IN_GENOME[0]}
BASE_GENOME=${FILE_GENOME#${INPUT_DIR}/}
BASE_GENOME=${BASE_GENOME%.fasta}

# query genome
INPUT_DIR=/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1
IN_QUERY=($(find $INPUT_DIR -maxdepth 1 -name "MesAur1.0.fa"))
FILE_QUERY=${IN_QUERY[0]}
BASE_QUERY=${FILE_QUERY#${INPUT_DIR}/}
BASE_QUERY=${BASE_QUERY%.fa}
BASE_QUERY="Siomi"

# query .gtf
INPUT_DIR=.
IN_GTF=($(find $INPUT_DIR -maxdepth 1 -name "*gff3"))
FILE_GTF=${IN_GTF[0]}
BASE_GTF=${FILE_GTF#${INPUT_DIR}/}
BASE_GTF=${BASE_GTF%.gff3}

# script
SCRIPT=/common/WORK/fhorvat/programi/python/packages/bin/liftoff

# ----------------Commands------------------- #
# extract introns and coverage from bam file
$SCRIPT \
-g ${FILE_GTF} \
-p ${THREADS} \
-o ./${BASE_GTF}.${BASE_QUERY}.liftoff.gff \
-f ./feature_list.txt \
${FILE_GENOME} \
${FILE_QUERY}
