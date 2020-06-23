#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=10g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -j oe
#PBS -N pbs.annotation_pwk
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# set names and date
NAME=mouse
SCI_NAME=mus_musculus_pwkphj
ENSEMBL_NAME=mpwkphj
ENSEMBL_RELEASE=93
GENOME_USCS=PWK_PhJ_v1
GENOME_ENSEMBL=PWK_PhJ_v1
ASSEMBLY_ID=GCA_001624775.1
ASSEMBLY_REPORT=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/624/775/GCA_001624775.1_PWK_PhJ_v1/GCA_001624775.1_PWK_PhJ_v1_assembly_report.txt

# create dir
DATE=`date +%Y%m%d`
OUTPATH=/common/WORK/fhorvat/annotation
OUTDIR=${OUTPATH}/${NAME}/${GENOME_USCS}.${GENOME_ENSEMBL}.${ASSEMBLY_ID}
mkdir -p $OUTDIR
cd $OUTDIR

## set URLs
# UCSC
UCSC_FTP=http://hgdownload.soe.ucsc.edu/hubs/mouseStrains
UCSC_URL=${UCSC_FTP}/${ASSEMBLY_ID}_${GENOME_ENSEMBL}

# ENSEMBL
ENSEMBL_URL=ftp://ftp.ensembl.org/pub/release-$ENSEMBL_RELEASE
ENSEMBL_GTF=`curl -l ${ENSEMBL_URL}/gtf/${SCI_NAME}/ | grep "gtf.gz" | grep -v "abinitio\|chr"`
#ENSEMBL_FASTA=`curl -l ${ENSEMBL_URL}/fasta/${SCI_NAME}/dna/ | grep "dna.toplevel.fa.gz"`

# get links
LINKS=(\
${UCSC_URL}/${ASSEMBLY_ID}_${GENOME_ENSEMBL}.2bit \
${UCSC_URL}/${ASSEMBLY_ID}_${GENOME_ENSEMBL}.rmsk.fa.out.gz \
${UCSC_URL}/${ASSEMBLY_ID}_${GENOME_ENSEMBL}.chrom.sizes
${ENSEMBL_URL}/gtf/${SCI_NAME}/${ENSEMBL_GTF} \
#${ENSEMBL_URL}/fasta/${SCI_NAME}/dna/${ENSEMBL_FASTA} \
${ASSEMBLY_REPORT} \
)

# ----------------Commands------------------- #
#### download
for LINK in "${LINKS[@]}"; do wget $LINK; done

#### clean
## repeatMasker
RMSK=(`ls *fa.out.gz`)
RMSK_NEW=rmsk.${GENOME_USCS}.${DATE}.raw.fa.out.gz
mv $RMSK $RMSK_NEW

## ensembl
mv ${ENSEMBL_GTF} ensembl.${ENSEMBL_RELEASE}.${GENOME_ENSEMBL}.${DATE}.gtf.gz
#mv ${ENSEMBL_FASTA} ensembl.${ENSEMBL_RELEASE}.${GENOME_ENSEMBL}.fa.gz

## genome
# convert 2bit to fasta, index the genome, gzip genome
mv ${ASSEMBLY_ID}_${GENOME_ENSEMBL}.chrom.sizes ${GENOME_USCS}.chrom.sizes
mv ${ASSEMBLY_ID}_${GENOME_ENSEMBL}.2bit ${GENOME_USCS}.2bit
twoBitToFa ${GENOME_USCS}.2bit ${GENOME_USCS}.fa
samtools faidx ${GENOME_USCS}.fa
gzip ${GENOME_USCS}.fa

## clean using custom R script
cp /common/WORK/fhorvat/annotation/scripts/02a.clean_files.R $OUTDIR
Rscript --vanilla 02a.clean_files.R $ENSEMBL_RELEASE $ENSEMBL_NAME $OUTDIR

## copy scripts to new dir
cp /common/WORK/fhorvat/annotation/scripts/01b.download_annotation.mouse_strains.pbs.sh $OUTDIR
