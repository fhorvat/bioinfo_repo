#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=10g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.03.trinotate_sqlite
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=1
MEMORY=10G

IN_DIR=/common/WORK/fhorvat/programi/Trinotate/Trinotate-Trinotate-v3.2.0
SCRIPT=$IN_DIR/Trinotate

# ----------------Commands------------------- #
# add to database
$SCRIPT $IN_DIR/admin/Trinotate.sqlite init \
--gene_trans_map s_GV_Mov10l_KO.DN.Trinity.l1000.fasta.gene_trans_map \
--transcript_fasta s_GV_Mov10l_KO.DN.Trinity.l1000.fasta \
--transdecoder_pep ./s_GV_Mov10l_KO.DN.Trinity.l1000.fasta.transdecoder_dir/longest_orfs.pep

# load BLAST homologies
$SCRIPT $IN_DIR/admin/Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6
$SCRIPT $IN_DIR/admin/Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6

# load Pfam domain entries
$SCRIPT $IN_DIR/admin/Trinotate.sqlite LOAD_pfam TrinotatePFAM.out

# output an Annotation Report
$SCRIPT $IN_DIR/admin/Trinotate.sqlite report > trinotate_annotation_report.xls
