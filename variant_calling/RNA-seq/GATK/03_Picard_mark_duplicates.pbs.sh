#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N Picard
#PBS -l select=ncpus=6:mem=40g
#PBS -J 0-4
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
INPUT_DIR=/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Analysis/variant_calling/GATK/02_STAR_2pass_map
IN_SEQ=(`ls ${INPUT_DIR}/*.sam`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE##*$INPUT_DIR/}
BASE=${BASE%_Aligned.out.sam}

PICARD_PATH=/common/WORK/fhorvat/programi/Picard/picard.jar

# ----------------Commands------------------- #
# add read groups, sort
java -jar $PICARD_PATH AddOrReplaceReadGroups I=$FILE O=${BASE}_sorted.bam SO=coordinate RGID=${BASE}_gatk RGLB=library RGPL=illumina RGPU=platform RGSM=$BASE

# mark duplicates, create index
java -jar $PICARD_PATH MarkDuplicates I=${BASE}_sorted.bam O=${BASE}_deduped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${BASE}.metrics.log
