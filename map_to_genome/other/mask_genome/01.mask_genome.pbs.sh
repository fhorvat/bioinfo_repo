#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.mask_genome
#PBS -l select=ncpus=10:mem=20g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# parameters
THREADS=10

# dirs and files
IN_FASTA=/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/mm10.fa.gz
IN_BED=/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/Documentation/L1s_nested_ours_20180516.ZJM.tidy.mm10_diff.bed

# ----------------Commands------------------- #
# ungzip
unpigz -p $THREADS $IN_FASTA

# mask fasta
bedtools maskfasta -fi ${IN_FASTA%.gz} -bed $IN_BED -fo mm10.L1s_nested_ours.diff.fa
samtools faidx mm10.L1s_nested_ours.diff.fa

# gzip
pigz -p $THREADS mm10.L1s_nested_ours.diff.fa
pigz -p $THREADS ${IN_FASTA%.gz}
