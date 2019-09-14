#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=8:mem=600g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.bridger_ts_assembly
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=8
MEMORY=590G

BRIDGER=/common/WORK/fhorvat/programi/Bridger/Bridger_r2014-12-01/Bridger.pl

BASE=`basename ${PWD/bridger/}`
IN_DIR=../../../Raw/$BASE/Cleaned
FILES=(`ls ${IN_DIR}/s_*all.PE_[1,2].txt.gz`)

# ----------------Commands------------------- #
# run script
mkdir bridger_out_dir
perl $BRIDGER --seqType fq --left ${FILES[0]} --right ${FILES[1]} --CPU $THREADS --debug &> bridger_${BASE}.log

