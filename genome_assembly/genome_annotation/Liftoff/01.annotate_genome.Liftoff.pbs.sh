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
IN_GENOME=($(find $INPUT_DIR -name "*.fasta"))
FILE_GENOME=${IN_GENOME[0]}
BASE_GENOME=${FILE_GENOME#${INPUT_DIR}/}
BASE_GENOME=${BASE_GENOME%.fasta}

# query genome
INPUT_DIR=/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed/annotation/NCBI_references
IN_QUERY=($(find -L $INPUT_DIR -name "*genomic.fna"))
FILE_QUERY=${IN_QUERY[0]}
FILE_GTF=${FILE_QUERY/.fna/.gff}
BASE_QUERY=${FILE_QUERY#${INPUT_DIR}/}
BASE_QUERY=${BASE_QUERY%_genomic.fna}

# script
SCRIPT=/common/WORK/fhorvat/programi/python/Python3/bin/liftoff

# ----------------Commands------------------- #
# extract introns and coverage from bam file
$SCRIPT \
-t ${FILE_GENOME} \
-r ${FILE_QUERY} \
-g ${FILE_GTF} \
-p ${THREADS} \
-o ./${BASE_GENOME}.${BASE_QUERY}.liftoff.gff
