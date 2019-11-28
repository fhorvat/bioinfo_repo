#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03.blastn.pbs.sh
#PBS -l select=ncpus=6:mem=30g
#PBS -J 0-1
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# number of threads
THREADS=6

# blastn subject 
IN_DB=Sirena1_pseudo.Elob.Elobl

# blastn query input
INPUT_DIR=.
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 -name "s_*.fa"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fa}

# ----------------Commands------------------- #
# blast
blastn -task blastn-short -query $FILE -db $IN_DB -out ${BASE}.sam -num_threads $THREADS -outfmt 17 \
-perc_identity 80 \
-word_size 4 \
-evalue 10 \
-dust no \
-soft_masking false \
-gapopen 10 \
-gapextend 20 \
-penalty -1
