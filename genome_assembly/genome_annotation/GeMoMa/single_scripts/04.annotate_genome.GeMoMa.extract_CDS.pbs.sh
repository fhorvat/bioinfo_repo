#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=10g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.04.annotate_genome.GeMoMa
#PBS -J 0-1
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1
MEMORY=10G

INPUT_DIR=./NCBI_references
IN_GENOME=($(find $INPUT_DIR -name "*genomic.fna"))
FILE_GENOME=${IN_GENOME[$PBS_ARRAY_INDEX]}
BASE_GENOME=${FILE_GENOME#${INPUT_DIR}/}
BASE_GENOME=${BASE_GENOME%.fna}

FILE_GTF=${FILE_GENOME/.fna/.gff}

SCRIPT=/common/WORK/fhorvat/programi/GeMoMa/GeMoMa-1.6.4/GeMoMa-1.6.4.jar

# ----------------Commands------------------- #
# create output dir
mkdir CDS.${BASE_GENOME}

# extract introns and coverage from bam file
java -jar $SCRIPT CLI Extractor \
a=$FILE_GTF \
g=$FILE_GENOME \
p=true \
c=true \
outdir=CDS.${BASE_GENOME}
