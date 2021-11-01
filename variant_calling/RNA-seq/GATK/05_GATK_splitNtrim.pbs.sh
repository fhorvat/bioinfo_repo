#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N splitNTrim
#PBS -l select=ncpus=6:mem=40g
#PBS -J 0-4
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
INPUT_DIR=/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Analysis/variant_calling/GATK/03_Picard
IN_SEQ=(`ls ${INPUT_DIR}/*.bam`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE##*$INPUT_DIR/}
BASE=${BASE%_deduped.bam}

GENOME_PATH=/common/WORK/fhorvat/reference/mouse/mm10/UCSC/mm10.fa
SCRIPT_PATH=/common/WORK/fhorvat/programi/GATK/GenomeAnalysisTK.jar

# ----------------Commands------------------- #
# split 
java -jar $SCRIPT_PATH -T SplitNCigarReads -R $GENOME_PATH -I $FILE -o ${BASE}_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

