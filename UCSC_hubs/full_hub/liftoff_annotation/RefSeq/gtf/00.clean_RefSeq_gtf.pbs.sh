#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=10g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.00.clean_RefSeq_gtf
#PBS -J 0-4
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
INPUT_DIR=.
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 -name "*_genomic.gtf"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}

# ----------------Commands------------------- #
# remove empty transcript_id fields 
sed -ie 's/transcript_id ""; //g' $FILE
