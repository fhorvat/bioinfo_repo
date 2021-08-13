#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.blastn.pbs.sh
#PBS -l select=ncpus=6:mem=30g
#PBS -J 0-1
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
IN_DB=${INPUT_DIR}/blastn_database/${BASE_GENOME}

# blastn query input
INPUT_DIR=.
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 -name "*.consensus.fasta"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fasta}

# ----------------Commands------------------- #
# blast
blastn -task rmblastn -query ${FILE} -db $IN_DB -out ${BASE}.blastn.txt -num_threads ${THREADS} -outfmt 6 \
-perc_identity 90 \
-dust no \
-soft_masking false

#-gapopen 10 \
#-gapextend 20 \
#-penalty -1
