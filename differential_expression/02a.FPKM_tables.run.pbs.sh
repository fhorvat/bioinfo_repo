#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.FPKM_tables
#PBS -l select=ncpus=1:mem=20g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# set variables
SCRIPT=02b.FPKM_tables.script.R
ENSEMBL_VERSION="93"
GENOME_PATH="/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"
MAPPED_DIR="STAR_mm10"
THREADS="1"
GROUPING_VARIABLES="stage genotype"

# experiment
EXPERIMENT=${PWD%/Analysis*}
EXPERIMENT=${EXPERIMENT##*/}

# mapped path
MAPPED_PATH=${PWD%/Analysis*}
MAPPED_PATH=${MAPPED_PATH}"/Data/Mapped/"${MAPPED_DIR}

# ----------------Commands------------------- #
# run script
Rscript $SCRIPT \
--experiment $EXPERIMENT \
--ensembl_version $ENSEMBL_VERSION \
--genome_path $GENOME_PATH \
--mapped_path $MAPPED_PATH \
--threads $THREADS \
--grouping_variables $GROUPING_VARIABLES
