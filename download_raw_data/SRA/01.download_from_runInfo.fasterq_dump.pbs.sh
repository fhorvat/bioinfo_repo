#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=4:mem=10g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.fasterq_dump
#PBS -j oe
#PBS -J 0-2
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=4

OUT_PATH=.
SAMPLE_TABLE=./*runInfo.txt
SRA_LIST=(`awk -F ',' -v col=Run 'NR==1{for(i=1;i<=NF;i++){if($i==col){c=i;break}} print $c} NR>1{print $c}' $SAMPLE_TABLE | sed '1d'`)
echo ${#SRA_LIST[@]}

SRA_RUN=${SRA_LIST[$PBS_ARRAY_INDEX]}

# ----------------Commands------------------- #
# download data
fasterq-dump \
--outdir $OUT_PATH \
--mem 9500MB \
--temp /common/WORK/fhorvat/tmp/sra \
--threads $THREADS \
--split-3 \
--skip-technical \
$SRA_RUN

# gzip files
for fq_file in `find $OUT_PATH -name "${SRA_RUN}*fastq"`; do pigz -p $THREADS $fq_file; done
