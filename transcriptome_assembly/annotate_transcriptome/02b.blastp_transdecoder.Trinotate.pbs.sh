#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02b.blastp_Trinotate
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

UNIPROT=/common/WORK/fhorvat/programi/Trinotate/Trinotate-Trinotate-v3.2.0/admin/uniprot_sprot.pep

# ----------------Commands------------------- #
# search Trinity transcripts against UniProt
blastp -query $FILE -db $UNIPROT -num_threads $THREADS -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastp.outfmt6
