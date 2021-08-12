#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=10g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.create_search_index
#PBS -J 0-1
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
INPUT_DIR=.
IN_SEQ=($(find $INPUT_DIR -name "*liftoff.gff"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.gff}

# ----------------Commands------------------- #
# required for indexing step
grep -v "^#" ${BASE}.info.txt | awk '{printf "%s\t%s,%s,%s,%s,%s\n", $1,$2,$3,$8,$9,$10}' > ${BASE}.nameIndex.txt

# create index
ixIxx ${BASE}.nameIndex.txt ${BASE}.ix ${BASE}.ixx

# set permissions
chmod 744 ${BASE}.ix ${BASE}.ixx
