#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=100g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01b.annotate_genome.GeMoMa
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1
MEMORY=100G

INPUT_DIR=/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.RNAseq/Data/Mapped/STAR_Siomi/5_merged_replicates
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 -name "*.bam" -and -name "*all*"))
FILE=${IN_SEQ[0]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}

SCRIPT=/common/WORK/fhorvat/programi/GeMoMa/GeMoMa-1.6.4/GeMoMa-1.6.4.jar

# ----------------Commands------------------- #
# create output dir
mkdir rna_seq.evidence.${BASE}

# extract introns and coverage from bam file
java -jar $SCRIPT CLI ERE \
m=${FILE} \
s=FR_FIRST_STRAND \
c=true \
outdir=rna_seq.evidence.${BASE}
