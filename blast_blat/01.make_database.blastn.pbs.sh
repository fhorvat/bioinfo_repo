#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.make_database.blastn.pbs.sh
#PBS -l select=ncpus=1:mem=1g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
INPUT_DIR=../
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 -name "*.fasta"))
FILE=${IN_SEQ[0]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fasta}

# ----------------Commands------------------- #
makeblastdb -in $IN_SEQ -out ./${BASE} -dbtype nucl
