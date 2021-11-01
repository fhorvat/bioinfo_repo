#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N GATK_HC
#PBS -l select=ncpus=6:mem=40g
#PBS -J 0-4
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=4

INPUT_DIR=/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Analysis/variant_calling/GATK/05_BQSR
IN_SEQ=(`ls ${INPUT_DIR}/*.bam`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE##*$INPUT_DIR/}
BASE=${BASE%_recal.bam}

GENOME_PATH=/common/WORK/fhorvat/reference/mouse/mm10/UCSC/mm10.fa
SCRIPT_PATH=/common/WORK/fhorvat/programi/GATK/GenomeAnalysisTK.jar

# ----------------Commands------------------- #
# variant calling
java -jar $SCRIPT_PATH -T HaplotypeCaller -R $GENOME_PATH -I $FILE -dontUseSoftClippedBases -stand_call_conf 20.0 -o ${BASE}.vcf -nct $THREADS
