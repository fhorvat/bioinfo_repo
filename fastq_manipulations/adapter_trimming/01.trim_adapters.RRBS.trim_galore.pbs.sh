#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=10g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.trim_adapters.RRBS
#PBS -j oe
#PBS -J 0-3
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=60g

IN_DIR=../Links
IN_SEQ=($(find $IN_DIR -maxdepth 1 -name "*.txt.gz" | grep -v "all"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.txt.gz}

# ----------------Commands------------------- #
# right trim
trim_galore --gzip --cores 4 --rrbs --basename ${BASE} ${FILE}
