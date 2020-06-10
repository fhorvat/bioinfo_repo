#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=250g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02b.ts_assembly.bridger
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=250G

BRIDGER=/common/WORK/fhorvat/programi/Bridger/Bridger_r2014-12-01/Bridger.pl

# ----------------Commands------------------- #
# run transcriptome assembly
NCPU=12
INFOLDER=.
OUTFOLDER=.
PREFIX=19_20_merged
perl $BRIDGER --seqType fq --left $INFOLDER/${PREFIX}_1.fq.gz --right $INFOLDER/${PREFIX}_2.fq.gz --CPU $NCPU --debug &> bridger_${PREFIX}.log

