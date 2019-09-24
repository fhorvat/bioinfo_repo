#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01a.summarizeOverlaps
#PBS -l select=ncpus=10:mem=150g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #

# ----------------Commands------------------- #
# run R script which summarizes overlaps between exons and reads in bam files
Rscript 01b.summarizeOverlaps.pbs.R
