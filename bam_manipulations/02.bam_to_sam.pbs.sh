#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.bam_to_sam
#PBS -l select=ncpus=12:mem=50g
#PBS -j oe
#PBS -J 0-2
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# number of threads
THREADS=12

# files
INPUT_DIR=/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.deduplicated.smallRNAseq/Data/Mapped/STAR_Siomi.multimappers
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 -name "*.bam" -and -name "*18to32nt*"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}

# ----------------Commands------------------- #
# bam to sam
samtools view -@ ${THREADS} -h -o ${BASE}.sam ${FILE}
