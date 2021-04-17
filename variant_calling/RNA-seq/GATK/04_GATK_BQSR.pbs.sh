#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N BQSR
#PBS -l select=ncpus=6:mem=40g
#PBS -J 0-4
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=4

INPUT_DIR=/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Analysis/variant_calling/GATK/04_SplitNTrim
IN_SEQ=(`ls ${INPUT_DIR}/*.bam`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE##*$INPUT_DIR/}
BASE=${BASE%_split.bam}

GENOME_PATH=/common/WORK/fhorvat/reference/mouse/mm10/UCSC/mm10.fa
SNP_VCF_PATH=/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Analysis/variant_calling/GATK/05_BQSR/mgp.v3.snps.rsIDdbSNPv137_UCSC_seqnames.vcf.gz
SCRIPT_PATH=/common/WORK/fhorvat/programi/GATK/GenomeAnalysisTK.jar

# ----------------Commands------------------- #
# create tables with before/after base quality recalibration
#java -jar $SCRIPT_PATH -T BaseRecalibrator -R $GENOME_PATH -I $FILE -knownSites $SNP_VCF_PATH -o ${BASE}_data.table -nct $THREADS
#java -jar $SCRIPT_PATH -T BaseRecalibrator -R $GENOME_PATH -I $FILE -knownSites $SNP_VCF_PATH -BQSR ${BASE}_data.table -o ${BASE}_recal_data.table -nct $THREADS

# optional - plot results
#java -jar $SCRIPT_PATH -T AnalyzeCovariates -R $GENOME_PATH -before ${BASE}_data.table -after ${BASE}_recal_data.table -plots ${BASE}_plots.pdf

# recalibrate base quality
java -jar $SCRIPT_PATH -T PrintReads -R $GENOME_PATH -I $FILE -BQSR ${BASE}_recal_data.table -o ${BASE}_recal.bam -nct $THREADS
