#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03.gene_trans_map
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1
MEMORY=20G

INPUT_DIR=.
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 -name "s_GV_Mov10l_KO.DN.Trinity.l1000.fasta"))
FILE=${IN_SEQ[0]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fasta}

SCRIPT=/common/WORK/fhorvat/programi/Trinity/trinityrnaseq-Trinity-v2.8.4/util/support_scripts/get_Trinity_gene_to_trans_map.pl

# ----------------Commands------------------- #
# create gene-transcript relationships
$SCRIPT $FILE > ${BASE}.fasta.gene_trans_map
