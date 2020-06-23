#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=10g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -j oe
#PBS -N pbs.annotation_cow
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# set names and date
NAME=mouse
SCI_NAME=mus_musculus
ENSEMBL_NAME=mmusculus
ENSEMBL_RELEASE=99
GENOME_USCS=mm10
GENOME_ENSEMBL=GRCm38.p6
ASSEMBLY_ID=GCA_000001635.2
ASSEMBLY_REPORT=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.20_GRCm38/GCF_000001635.20_GRCm38_assembly_report.txt

# date and outdir
DATE=`date +%Y%m%d`
OUTPATH=/common/DB/genome_reference
OUTDIR=${OUTPATH}/${NAME}/${GENOME_USCS}.${GENOME_ENSEMBL/.p6/}.${ASSEMBLY_ID}
OUTDIR=$PWD

# ENSEMBL
ENSEMBL_URL=ftp://ftp.ensembl.org/pub/release-$ENSEMBL_RELEASE
ENSEMBL_GTF=`curl -l ${ENSEMBL_URL}/gtf/${SCI_NAME}/ | grep "gtf.gz" | grep -v "abinitio\|chr"`

# get links
LINKS=(\
${ENSEMBL_URL}/gtf/${SCI_NAME}/${ENSEMBL_GTF} \
)

# ----------------Commands------------------- #
#### download
for LINK in "${LINKS[@]}"; do wget $LINK; done

## ensembl
mv ${ENSEMBL_GTF} ensembl.${ENSEMBL_RELEASE}.${GENOME_ENSEMBL}.${DATE}.gtf.gz

## clean using ustom R script
Rscript 03b.clean_files.ensembl_version.R $ENSEMBL_RELEASE $ENSEMBL_NAME $OUTDIR
