#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02a.blastx_Trinotate
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6
MEMORY=20G

INPUT_DIR=.
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 -name "*l1000.fasta"))
FILE=${IN_SEQ[0]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.fasta}

UNIPROT=/common/WORK/fhorvat/programi/Trinotate/Trinotate-Trinotate-v3.2.0/admin/uniprot_sprot.pep

# ----------------Commands------------------- #
# search Trinity transcripts against UniProt
blastx -query $FILE -db $UNIPROT -num_threads $THREADS -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastx.outfmt6
