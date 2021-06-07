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
NAME=rat
SCI_NAME=rattus_norvegicus
ENSEMBL_NAME=bnorvegicus
ENSEMBL_RELEASE=99
GENOME_USCS=rn7
GENOME_ENSEMBL=rn7
ASSEMBLY_ID=GCA_015227675.2
ASSEMBLY_REPORT=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/227/675/GCF_015227675.2_mRatBN7.2/GCF_015227675.2_mRatBN7.2_assembly_report.txt

# create dir
DATE=`date +%Y%m%d`

## set URLs
# UCSC
UCSC_FTP=ftp://hgdownload.soe.ucsc.edu/goldenPath
UCSC_URL=${UCSC_FTP}/${GENOME_USCS}

# RefSeq
REFSEQ_GFF_LINK=${ASSEMBLY_REPORT/_assembly_report.txt/_genomic.gff.gz}

# get links
LINKS=(\
${UCSC_URL}/bigZips/${GENOME_USCS}.2bit \
${UCSC_URL}/bigZips/${GENOME_USCS}.fa.out.gz \
${UCSC_URL}/bigZips/${GENOME_USCS}.chrom.sizes
${ASSEMBLY_REPORT} \
${REFSEQ_GFF_LINK} \
)

# ----------------Commands------------------- #
#### download
for LINK in "${LINKS[@]}"; do wget $LINK; done

#### clean
## repeatMasker
RMSK=($(find . -name "*fa.out.gz"))
RMSK_NEW=rmsk.${GENOME_USCS}.${DATE}.raw.fa.out.gz
mv $RMSK $RMSK_NEW

## genome
# convert 2bit to fasta, index the genome, bgzip genome
twoBitToFa ${GENOME_USCS}.2bit ${GENOME_USCS}.fa
bgzip ${GENOME_USCS}.fa
samtools faidx ${GENOME_USCS}.fa.gz
