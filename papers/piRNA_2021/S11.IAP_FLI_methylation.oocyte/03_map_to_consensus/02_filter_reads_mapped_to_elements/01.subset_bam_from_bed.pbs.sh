#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.filter_bam
#PBS -l select=ncpus=1:mem=30g
#PBS -j oe
#PBS -J 0-2
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# files
INPUT_DIR=/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.bisulfite/Data/Mapped/Bismark_Siomi.trimmed
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 -name "*.deduplicated.bam"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}

BED=/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/IAP/annotation/IAP.potentially_young.ordered_by_ORFs.20201031.IAPLTR3.bed

# ----------------Commands------------------- #
# filter by bed
samtools view -b -L $BED $FILE > $BASE.bam

# index
samtools index ${BASE}.bam
