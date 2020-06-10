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
INPUT_DIR=.
IN_SEQ=${INPUT_DIR}/Sirena1_pseudo.Elob.Elobl.fa
BASE=${IN_SEQ#${INPUT_DIR}/}
BASE=${BASE%.fa}

# ----------------Commands------------------- #
makeblastdb -in $IN_SEQ -out $BASE -dbtype nucl -parse_seqids

