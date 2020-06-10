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
NAME=cow
SCI_NAME=bos_taurus
ENSEMBL_NAME=btaurus
ENSEMBL_RELEASE=96
GENOME_USCS=bosTau9
GENOME_ENSEMBL=ARS-UCD1.2
ASSEMBLY_ID=GCA_002263795.2
ASSEMBLY_REPORT=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/263/795/GCF_002263795.1_ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_assembly_report.txt

# create dir
DATE=`date +%Y%m%d`
OUTPATH=/common/DB/genome_reference
OUTDIR=${OUTPATH}/${NAME}/${GENOME_USCS}.${GENOME_ENSEMBL}.${ASSEMBLY_ID}
mkdir -p $OUTDIR
cd $OUTDIR

## set URLs
# UCSC
UCSC_FTP=ftp://hgdownload.soe.ucsc.edu/goldenPath
UCSC_URL=${UCSC_FTP}/${GENOME_USCS}

# ENSEMBL
ENSEMBL_URL=ftp://ftp.ensembl.org/pub/release-$ENSEMBL_RELEASE
ENSEMBL_GTF=`curl -l ${ENSEMBL_URL}/gtf/${SCI_NAME}/ | grep "gtf.gz" | grep -v "abinitio\|chr"`
#ENSEMBL_FASTA=`curl -l ${ENSEMBL_URL}/fasta/${SCI_NAME}/dna/ | grep "dna.toplevel.fa.gz"`

# get links
LINKS=(\
${UCSC_URL}/bigZips/${GENOME_USCS}.2bit \
${UCSC_URL}/bigZips/${GENOME_USCS}.fa.out.gz \
${UCSC_URL}/bigZips/${GENOME_USCS}.chrom.sizes
#${UCSC_URL}/database/refGene.txt.gz \
${ENSEMBL_URL}/gtf/${SCI_NAME}/${ENSEMBL_GTF} \
#${ENSEMBL_URL}/fasta/${SCI_NAME}/dna/${ENSEMBL_FASTA} \
${ASSEMBLY_REPORT} \
)

# ----------------Commands------------------- #
#### download
for LINK in "${LINKS[@]}"; do wget $LINK; done

#### clean
## repeatMasker
RMSK=($(find . -name "*fa.out.gz"))
RMSK_NEW=rmsk.${GENOME_USCS}.${DATE}.raw.fa.out.gz
mv $RMSK $RMSK_NEW

## refGene - rename, convert to genePred to gtf, remove genePred
#mv refGene.txt.gz refGene.${GENOME_USCS}.${DATE}.txt.gz
#zcat refGene.${GENOME_USCS}.${DATE}.txt.gz | cut -f 2- | genePredToGtf file stdin stdout -source=refSeq.${DATE} | sort -k1,1 -k4,4 | gzip > refGene.${GENOME_USCS}.${DATE}.gtf.gz
#[ -f "refGene.${GENOME_USCS}.${DATE}.gtf.gz" ] && rm refGene.${GENOME_USCS}.${DATE}.txt.gz

## ensembl
mv ${ENSEMBL_GTF} ensembl.${ENSEMBL_RELEASE}.${GENOME_ENSEMBL}.${DATE}.gtf.gz
#mv ${ENSEMBL_FASTA} ensembl.${ENSEMBL_RELEASE}.${GENOME_ENSEMBL}.fa.gz

## genome
# convert 2bit to fasta, index the genome, bgzip genome
twoBitToFa ${GENOME_USCS}.2bit ${GENOME_USCS}.fa
bgzip ${GENOME_USCS}.fa
samtools faidx ${GENOME_USCS}.fa.gz

## clean using custom R script
cp /common/WORK/fhorvat/annotation/scripts/02a.clean_files.R $OUTDIR
Rscript --vanilla 02a.clean_files.R $ENSEMBL_RELEASE $ENSEMBL_NAME $OUTDIR

## copy scripts to new dir
cp /common/WORK/fhorvat/annotation/scripts/01a.download_annotation.common.pbs.sh $OUTDIR
