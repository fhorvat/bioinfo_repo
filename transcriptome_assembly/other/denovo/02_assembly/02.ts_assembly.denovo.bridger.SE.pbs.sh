#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=250g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.ts_assembly.bridger
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=250G

INPUT_DIR=..
BASE="pwd"

SCRIPT=/common/WORK/fhorvat/programi/Bridger/Bridger_r2014-12-01/Bridger.pl
SCRIPT_PAR="--seqType fq \
--clean \
--CPU $THREADS"

# ----------------Commands------------------- #
# concatenate fastq files to one file
cat $INPUT_DIR/*txt.gz > ${BASE}.merged.txt.gz

# create working dirs
mkdir -p bridger_out_dir/RawGraphs

# run transcriptome assembly
$SCRIPT $SCRIPT_PAR --single ${BASE}.merged.txt.gz &> ${BASE}.Bridger.log

# move .fasta from output dir
mv bridger_out_dir/Bridger.fasta .
mv Bridger.fasta ${BASE}.Bridger.fasta

# remove output dir
[ -f "${BASE}.Bridger.fasta" ] && rm -r bridger_out_dir

# get assembly statistics
/common/WORK/fhorvat/programi/Trinity/trinityrnaseq-Trinity-v2.5.1/util/TrinityStats.pl ${BASE}.Bridger.fasta > ${BASE}.Bridger.stats
 
