### INFO: R Script
### DATE: 24.03.2017. 
### AUTHOR: Filip Horvat
### PATH: 
rm(list = ls()); gc()
# options(bitmapType = 'cairo')

################################################################################### WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/1cell_motif_discovery")

################################################################################### LIBRARIES
library(data.table)
library(stringr)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(reshape2)

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(BiocParallel)
library(DESeq2)

################################################################################### PATH VARIABLES
inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mm10"

outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/1cell_motif_discovery/documentation"

mm10_gtf_path <- "/common/WORK/fhorvat/reference/mouse/mm10/EnsemblGenes/Ensembl_GRCm38.86.20161128.gtf.gz"

lib_path <- "/common/WORK/fhorvat/R_library"

################################################################################### SOURCE FILES
source(file.path(lib_path, "GffToGRanges.R"))

################################################################################### FUNCTIONS

################################################################################### SCRIPT PARAMS

################################################################################### TABLES
# experiment table
sample_table <- 
  read_delim("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/CNOT6L_library_size.txt", delim = "\t", col_names = T) %>% 
  dplyr::filter(str_detect(name, "1C")) %>% 
  dplyr::mutate(treatment = str_replace_all(name, "s_1C_|_X.*", ""))

# read in the ENSEMBL gene table, take protein coding genes
ensembl_gtf <- 
  read_delim(file = mm10_gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c")) %>% 
  dplyr::filter(str_detect(X9, "protein_coding"))

################################################################################### MAIN CODE
### take transcript with most exons for each gene
# get 5'UTR-s for each transcript, take only those longer than 7
gtf_5utr <- GffToGRanges(ensembl_gtf, "five_prime_utr")
gtf_5utr <- split(gtf_5utr, gtf_5utr$transcript_id)
gtf_5utr <- gtf_5utr[sum(width(gtf_5utr)) > 7]
gtf_5utr <- unlist(gtf_5utr)

# get exons for each transcript
gtf_trans <- GffToGRanges(ensembl_gtf, "exon")

# get unique gene_id/transcript_id combinations
gids <- unique(values(gtf_trans)[c('gene_id', 'transcript_id')])

# splits exon ranges on transcripts
gtf_trans <- split(gtf_trans, gtf_trans$transcript_id)

# orders transcripts based on number of exons in transcripts
gtf_trans <- gtf_trans[order(elementNROWS(gtf_trans), decreasing = T)] 

# keeps only first transcript of each gene (the one with most exons)
gtf_trans <- gtf_trans[!duplicated(gids$gene_id[match(names(gtf_trans), gids$transcript_id)])]

# take only transcripts which have 5' UTR
gtf_trans <- gtf_trans[names(gtf_trans) %in% gtf_5utr$transcript_id]

################################################################################### 
### get count of reads 
# counting overlaps
bamfiles <- Rsamtools::BamFileList(sample_table$sample_path, yieldSize = 2000000)
register(MulticoreParam())
se <- GenomicAlignments::summarizeOverlaps(features = gtf_trans,
                                           reads = bamfiles,
                                           mode = "Union",
                                           singleEnd = FALSE,
                                           ignore.strand = TRUE)

# set column data
colData(se) <- DataFrame(sample_table)

# DESeq2 differential expression - KO vs. WT
dds <- DESeqDataSet(se, design = ~treatment)
dds <- dds[rowSums(counts(dds)) > 1, ]
colnames(dds) <- sample_table$ID
dds <- DESeq(dds)

# get results 
res <- 
  results(dds, contrast = c("treatment", "KO", "WT")) %>% 
  as.data.frame(.) %>% 
  tibble::rownames_to_column(var = "ensemblID") 

################################################################################### 
# get IDs of top 100 transcripts by expression in 1-cell stage
top_transcriptIDs <- 
  res %>% 
  dplyr::top_n(n = 100, wt = log2FoldChange) %>% 
  dplyr::arrange(desc(log2FoldChange)) %$% 
  ensemblID %>% 
  as.character(.) 
  
# get ranges of 5'UTRs of top 100 transcripts, write ranges
top_transcript_5utr <-
  gtf_5utr[gtf_5utr$transcript_id %in% top_transcriptIDs] %>%
  as.data.frame(., row.names = NULL) %>%
  write_delim(path = file.path(outpath, "top100_KOvsWT_upreg_ensembl_genes_5UTR_coordinates.txt"))
