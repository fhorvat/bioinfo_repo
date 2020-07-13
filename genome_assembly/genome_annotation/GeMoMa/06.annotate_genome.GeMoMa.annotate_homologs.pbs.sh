#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=60g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.06.annotate_genome.GeMoMa
#PBS -J 0-1
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=20G

# tblastn results
INPUT_DIR=.
IN_TBLASTN=($(find $INPUT_DIR -name "tblastn_results.sorted.txt"))
FILE_TBLASTN=${IN_TBLASTN[$PBS_ARRAY_INDEX]}
BASE_TBLASTN=${FILE_TBLASTN#${INPUT_DIR}/}
BASE_TBLASTN=${BASE_TBLASTN%/*}

# CDS which were blasted to genomes
FILE_CDS=${INPUT_DIR}/"CDS".${BASE_TBLASTN}/cds-parts.fasta
FILE_ASSIGNMENT=${INPUT_DIR}/"CDS".${BASE_TBLASTN}/assignment.tabular
FILE_PROTEIN=${INPUT_DIR}/"CDS".${BASE_TBLASTN}/proteins.fasta

# RNAseq files - introns and coverage
FILE_INTRONS_1=${INPUT_DIR}/rna_seq.evidence.s_GoldHam_all/denoised_introns.gff
FILE_COVERAGE_1=${INPUT_DIR}/rna_seq.evidence.s_GoldHam_all/coverage.bedgraph

FILE_INTRONS_2=${INPUT_DIR}/rna_seq.evidence.s_testis_Mov10l_all/denoised_introns.gff
FILE_COVERAGE_2_FORWARD=${INPUT_DIR}/rna_seq.evidence.s_testis_Mov10l_all/coverage_forward.bedgraph
FILE_COVERAGE_2_REVERSE=${INPUT_DIR}/rna_seq.evidence.s_testis_Mov10l_all/coverage_reverse.bedgraph

# target genome
INPUT_DIR=../../
IN_GENOME=($(find $INPUT_DIR -name "*.fasta"))
FILE_GENOME=${IN_GENOME[0]}
BASE_GENOME=${FILE_GENOME#${INPUT_DIR}/}
BASE_GENOME=${BASE_GENOME%/*}

SCRIPT=/common/WORK/fhorvat/programi/GeMoMa/GeMoMa-1.6.4/GeMoMa-1.6.4.jar

# ----------------Commands------------------- #
# create output dir
mkdir gene_models.${BASE_TBLASTN}

# find homologous candidate transcripts in the target organism
java -jar $SCRIPT CLI GeMoMa \
s=${FILE_TBLASTN} \
t=${FILE_GENOME} \
c=${FILE_CDS} \
a=${FILE_ASSIGNMENT} \
q=${FILE_PROTEIN} \
i=${FILE_INTRONS_1} \
coverage=UNSTRANDED \
coverage_unstranded=${FILE_COVERAGE_1} \
i=${FILE_INTRONS_2} \
coverage=STRANDED \
coverage_forward=${FILE_COVERAGE_2_FORWARD} \
coverage_reverse=${FILE_COVERAGE_2_REVERSE} \
outdir=gene_models.${BASE_TBLASTN}
