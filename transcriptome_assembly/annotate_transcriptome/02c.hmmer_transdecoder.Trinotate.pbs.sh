#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02c.hmmscan_Trinotate
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6
MEMORY=20G

INPUT_DIR=./s_GV_Mov10l_KO.DN.Trinity.l1000.fasta.transdecoder_dir
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 -name "*pep"))
FILE=${IN_SEQ[0]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.pep}

PFAM=/common/WORK/fhorvat/programi/Trinotate/Trinotate-Trinotate-v3.2.0/admin/Pfam-A.hmm

# ----------------Commands------------------- #
# search through pfam database
hmmscan --cpu 12 --domtblout TrinotatePFAM.out $PFAM $FILE > pfam.log
