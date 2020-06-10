#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=1g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.fastqdump
#PBS -j oe
#PBS -J 0-8
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
OUT_PATH=.
SAMPLE_TABLE=./*runInfo.txt
SRA_LIST=(`awk -F '\t' -v col=Run 'NR==1{for(i=1;i<=NF;i++){if($i==col){c=i;break}} print $c} NR>1{print $c}' $SAMPLE_TABLE | sed '1d'`)
SRA_RUN=${SRA_LIST[$PBS_ARRAY_INDEX]}

# ----------------Commands------------------- #
fastq-dump --split-3 --gzip --origfmt --outdir $OUT_PATH $SRA_RUN
