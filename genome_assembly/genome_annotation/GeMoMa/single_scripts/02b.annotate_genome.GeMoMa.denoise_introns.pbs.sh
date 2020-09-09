#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02b.annotate_genome.GeMoMa
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1
MEMORY=20G

INPUT_DIR=./rna_seq.evidence.s_testis_Mov10l_all
IN_GTF=($(find $INPUT_DIR -name "introns.gff"))
FILE_GTF=${IN_GTF[0]}
BASE_GTF=${FILE_GTF#${INPUT_DIR}/}
BASE_GTF=${BASE_GTF%/*}

FILE_COVERAGE_FORWARD=${INPUT_DIR}/"coverage_forward.bedgraph"
FILE_COVERAGE_REVERSE=${INPUT_DIR}/"coverage_reverse.bedgraph"

SCRIPT=/common/WORK/fhorvat/programi/GeMoMa/GeMoMa-1.6.4/GeMoMa-1.6.4.jar

# ----------------Commands------------------- #
# extract introns and coverage from bam file
java -jar $SCRIPT CLI DenoiseIntrons \
i=$FILE_GTF \
c=STRANDED \
coverage_forward=$FILE_COVERAGE_FORWARD \
coverage_reverse=$FILE_COVERAGE_REVERSE \
outdir=${BASE_GTF}
