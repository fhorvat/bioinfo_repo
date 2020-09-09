#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=24:mem=80g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.annotate_genome.GeMoMa
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=24
MEMORY=80G

# target genome
INPUT_DIR=/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed
IN_GENOME=($(find $INPUT_DIR -name "*.fasta"))
FILE_GENOME=${IN_GENOME[0]}
BASE_GENOME=${FILE_GENOME#${INPUT_DIR}/}
BASE_GENOME=${BASE_GENOME%.fasta}

## query genomes
INPUT_DIR=./NCBI_references
IN_QUERY=($(find -L $INPUT_DIR -name "*genomic.fna"))
FILE_QUERY_1=${IN_QUERY[0]}
FILE_QUERY_2=${IN_QUERY[1]}
FILE_GTF_1=${FILE_QUERY_1/.fna/.gff}
FILE_GTF_2=${FILE_QUERY_2/.fna/.gff}

# RNA-seq files
INPUT_DIR=.
FILE_INTRONS_1=${INPUT_DIR}/rna_seq.evidence.s_GoldHam_all/introns.gff
FILE_INTRONS_2=${INPUT_DIR}/rna_seq.evidence.s_testis_Mov10l_all/introns.gff
FILE_COVERAGE_1=${INPUT_DIR}/rna_seq.evidence.s_GoldHam_all/coverage.bedgraph
FILE_COVERAGE_2_FORWARD=${INPUT_DIR}/rna_seq.evidence.s_testis_Mov10l_all/coverage_forward.bedgraph
FILE_COVERAGE_2_REVERSE=${INPUT_DIR}/rna_seq.evidence.s_testis_Mov10l_all/coverage_reverse.bedgraph

# script
SCRIPT=/common/WORK/fhorvat/programi/GeMoMa/GeMoMa-1.6.4/GeMoMa-1.6.4.jar

# ----------------Commands------------------- #
# extract introns and coverage from bam file
java -jar $SCRIPT CLI GeMoMaPipeline \
t=${FILE_GENOME} \
s="own" i="golden_hamster" a=${FILE_GTF_1} g=${FILE_QUERY_1} \
s="own" i="chinese_hamster" a=${FILE_GTF_2} g=${FILE_QUERY_2} \
tblastn=true \
r="EXTRACTED" \
introns=${FILE_INTRONS_1} coverage="UNSTRANDED" coverage_unstranded=${FILE_COVERAGE_1} \
AnnotationFinalizer.r="NO" \
pc=true \
outdir=./GeMoMa_results \
threads=${THREADS}
