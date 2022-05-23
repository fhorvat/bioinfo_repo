#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.translate_genome.transeq
#PBS -l select=ncpus=1:mem=1g
#PBS -j oe
 
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# number of threads
THREADS=1

# genome file
INPUT_DIR=.
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 -name "AriVul.fix.Dicer1.mRNA_genomic.fa"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fa}
BASE=${BASE%.fasta}

# script path
SCRIPT=/common/WORK/fhorvat/programi/emboss/EMBOSS-6.6.0/bin/transeq

# ----------------Commands------------------- #
# translate
${SCRIPT} -sequence ${FILE} -outseq ${BASE}.translated.fasta -frame=6 -clean
