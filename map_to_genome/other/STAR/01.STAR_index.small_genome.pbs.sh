#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=4:mem=5g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -j oe
#PBS -N pbs.01.STAR_index
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script parameters
THREADS=4

# in and out
IN_PATH=.
IN_GENOME=($(find $IN_PATH -name "*fasta"))
OUTDIR=.

# ----------------Commands------------------- #
# generate index (--genomeSAindexNbases = min(14, log2(GenomeLength)/2 - 1))
STAR --runThreadN $THREADS --runMode genomeGenerate --genomeDir $OUTDIR --genomeFastaFiles $IN_GENOME --genomeSAindexNbases 4
