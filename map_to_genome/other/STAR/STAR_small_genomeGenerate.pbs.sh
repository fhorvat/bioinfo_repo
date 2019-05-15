#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=10:mem=50g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N STAR
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
IN_GENOME=/common/WORK/fhorvat/reference/mouse/mm10/UCSC/rDNA_45s_mm10_BK000964.3.fasta
OUTDIR=/common/WORK/fhorvat/reference/mouse/mm10/genome_indexes/STAR/rDNA_45S_BK000964.3
THREADS=8

# ----------------Commands------------------- #
# generate index (--genomeSAindexNbases = min(14, log2(GenomeLength)/2 - 1))
STAR --runThreadN $THREADS --runMode genomeGenerate --genomeDir $OUTDIR --genomeFastaFiles $IN_GENOME --genomeSAindexNbases 4

