#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.07.create_clade_presence_absence_matrix.stk
#PBS -l select=ncpus=1:mem=1g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# input 
INPUT_DIR=.
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 \( -name "*.tre" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE}
BASE="RAxML_fastTreeSH_Support.Dicer1.Metazoa.OrthoDB.20211117.taxonomy_and_sequences_PS_CDD.ALL_msa.shsupport"

SCRIPT="/common/WORK/fhorvat/programi/Perl/libraries/bin/stk_create_matrix.pl"

# ----------------Commands------------------- #
# 
${SCRIPT} --dir ${INPUT_DIR} --output ${BASE}.nex --format Nexus
