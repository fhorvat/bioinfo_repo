#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N MOSAIK_CS
#PBS -J 0-7
#PBS -l select=ncpus=6:mem=50g
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
INPUT_DIR=/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/Cao_2014_pig/Data/Raw/Links
IN_SEQ=($INPUT_DIR/*_F3.csfasta.gz)
FILE_SEQ=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE_SEQ##*$INPUT_DIR/}
BASE=${BASE%_F3.csfasta.gz}
FILE_QUAL=${INPUT_DIR}/${BASE}_F3_QV.qual.gz

# ----------------Commands------------------- #
# build read binaries
MosaikBuild -fr $FILE_SEQ -fq $FILE_QUAL -out ${BASE}_cs.dat -st solid  &> ${BASE}_read_binaries.log 
