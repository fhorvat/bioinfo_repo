#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.quantify.Salmon
#PBS -l select=ncpus=8:mem=30g
#PBS -J 0-8
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=8
SALMON_INDEX=/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.RNAseq/Analysis/transcriptome_expression_ensembl/ensembl.93.MesAur1.0.20180920.cdna.all.salmon_index

INPUT_DIR=../../Data/Raw/Links
IN_SEQ=($(find ${INPUT_DIR} -name "*.txt.gz"))
UNIQ_SEQ=(`printf "%s\n" "${IN_SEQ[@]%.txt.gz}" | sort -u`)
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}

# ----------------Commands------------------- #
# Salmon quantify
salmon quant \
-i ${SALMON_INDEX} \
-l A \
-r ${FILE}.txt.gz \
-p ${THREADS} \
--validateMappings \
-o ${BASE}_quant
