#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.translate_genome.ORFs.IPS_translate
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
SCRIPT=/common/WORK/fhorvat/programi/InterProScan/interproscan-5.55-88.0/bin/nucleotide/translate

# ----------------Commands------------------- #
# translate
${SCRIPT} -i ${FILE} -o ${BASE}.IPS.translated.fasta
