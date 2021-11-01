#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=24:mem=80g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.annotate_genome.BRAKER
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=24
MEMORY=80G

IN_DIR_GEN=/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed/soft_masked_genome
IN_SEQ_GEN=($(find $IN_DIR_GEN -maxdepth 1 -name "*.fasta"))
FILE_GEN=${IN_SEQ_GEN[0]}

IN_DIR_RNASEQ_1=/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/Mapped/golden_hamster.Siomi
IN_DIR_RNASEQ_2=/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.RNAseq/Data/Mapped/STAR_Siomi
IN_SEQ_RNASEQ=$(find $IN_DIR_RNASEQ_1 $IN_DIR_RNASEQ_2 -maxdepth 1 -name "*.bam" |  paste -sd, -)
FILE_RNASEQ=${IN_SEQ_RNASEQ}

IN_DIR_PROT=/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed/annotation/protein_sequences
IN_SEQ_PROT=($(find $IN_DIR_PROT -maxdepth 1 -name "*.fasta"))
FILE_PROT=${IN_SEQ_PROT[0]}

BASE=${FILE_GEN#${IN_DIR_GEN}/}
BASE=${BASE%.fasta}

SCRIPT=/common/WORK/fhorvat/programi/BRAKER/BRAKER/scripts/braker.pl

# ----------------Commands------------------- #
# annotate the genome using RNA-seq and protein data
$SCRIPT \
--species=golden_hamster \
--genome=${FILE_GEN} \
--prot_seq=${FILE_PROT} \
--bam=${FILE_RNASEQ} \
--etpmode \
--UTR=off \
--softmasking \
--verbosity=4 \
--cores=${THREADS}
