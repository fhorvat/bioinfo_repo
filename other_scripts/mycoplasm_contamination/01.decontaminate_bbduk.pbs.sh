#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=30g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.decontaminate.bbduk
#PBS -j oe
#PBS -J 0-1
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6

IN_DIR=.
IN_SEQ=(`ls $IN_DIR/*.fastq.gz`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.fastq.gz}

REF=/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Data/Raw/mycoplasma_contamination/genomes/mycoplasma_genomes.fa

# ----------------Commands------------------- #
bbduk.sh in=$FILE k=18 ref=$REF threads=$THREADS 2> ${BASE}.k18.log
bbduk.sh in=$FILE k=21 ref=$REF threads=$THREADS 2> ${BASE}.k21.log
bbduk.sh in=$FILE k=25 ref=$REF threads=$THREADS 2> ${BASE}.k25.log
bbduk.sh in=$FILE k=27 ref=$REF threads=$THREADS 2> ${BASE}.k27.log
bbduk.sh in=$FILE k=31 ref=$REF threads=$THREADS 2> ${BASE}.k31.log

