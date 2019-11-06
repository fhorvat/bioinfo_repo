#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=2:mem=10g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N bbduk_adapter
#PBS -J 0-12
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
LEFT_ADAPTER=/common/WORK/fhorvat/Projekti/Svoboda/SMRT_oocytes_2017/Data/documentation/used_adapters.fasta
RIGHT_ADAPTER=/common/WORK/fhorvat/Projekti/Svoboda/SMRT_oocytes_2017/Data/documentation/used_adapters_revcomp.fasta
LEFT_BARCODES=/common/WORK/fhorvat/Projekti/Svoboda/SMRT_oocytes_2017/Data/documentation/used_barcodes.fasta
RIGHT_BARCODES=/common/WORK/fhorvat/Projekti/Svoboda/SMRT_oocytes_2017/Data/documentation/used_barcodes_revcomp.fasta

INPUT_DIR=/common/WORK/fhorvat/Projekti/Svoboda/SMRT_oocytes_2017/Data/Raw/Cleaned/02_demultiplexed/bbtools_seal
IN_SEQ=($INPUT_DIR/*.fastq)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE##*$INPUT_DIR/}
BASE=${BASE%%.fastq}

# ----------------Commands------------------- #
bbduk.sh in=$FILE out=${BASE}_1.fastq ref=$LEFT_ADAPTER rcomp=false ktrim=l k=13 mink=8 hdist=2 restrictleft=30
bbduk.sh in=${BASE}_1.fastq out=${BASE}_2.fastq ref=$RIGHT_ADAPTER rcomp=false ktrim=r k=12 mink=8 hdist=2 restrictright=30
bbduk.sh in=${BASE}_2.fastq out=${BASE}_3.fastq ref=$LEFT_BARCODES rcomp=false ktrim=l k=8 mink=4 hdist=2 restrictleft=16
bbduk.sh in=${BASE}_3.fastq out=${BASE}.fastq ref=$RIGHT_BARCODES rcomp=false ktrim=r k=7 mink=4 hdist=2 restrictright=16
rm ${BASE}_1.fastq ${BASE}_2.fastq ${BASE}_3.fastq
