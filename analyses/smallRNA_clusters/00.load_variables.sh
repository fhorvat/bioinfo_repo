#!/bin/bash
# ----------------Loading variables------------------- #
# set experiment
EXPERIMENT=DicerX_embryos

# base path
BASE_PATH=/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/DicerX_viral_infection/datasets/2019_Dec

# mapped and documentation
MAPPED_PATH=${BASE_PATH}/Data/Mapped/STAR_mm10/2_original_mapping
DOCUMENTATION_PATH=${BASE_PATH}/Data/Documentation

# genome and features
GENOME_PATH=/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2
INFO_PATH=${GENOME_PATH}/ensembl.93.GRCm38.p6.20180919.UCSCseqnames.geneInfo.csv
EXONS_PATH=${GENOME_PATH}/ensembl.93.GRCm38.p6.20180919.UCSCseqnames.reducedExons.RDS
RMSK_PATH=${GENOME_PATH}/rmsk.mm10.20180919.clean.fa.out.gz
MIRBASE_PATH=${GENOME_PATH}/miRBase.22.mm10.20181605.gff3

# other
THREADS=1
