#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.04a.blat_mRNA
#PBS -l select=ncpus=1:mem=10g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# number of threads
THREADS=6

# genome file
INPUT_DIR=/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed
IN_GENOME=($(find $INPUT_DIR -maxdepth 1 -name "*.fasta"))
FILE_GENOME=${IN_GENOME[0]}
BASE_GENOME=${FILE_GENOME#${INPUT_DIR}/}
BASE_GENOME=${BASE_GENOME%.fasta}

# blastn query input
INPUT_DIR=.
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 -name "*.mRNA.fasta"))
FILE=${IN_SEQ[0]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fasta}

# ----------------Commands------------------- #
# blat
blat -stepSize=5 -repMatch=1024 -minScore=0 -minIdentity=0 ${FILE_GENOME} ${FILE} ${BASE}.psl

# get score 
pslScore ${BASE}.psl > ${BASE}.score.psl
