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
sample_table <- read_delim("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/CNOT6L_library_size.txt", delim = "\t", col_names = T)

# get path to samples from table
sample_path <- 
  sample_table %>% 
  dplyr::filter(str_detect(name, "1C_WT")) %$%
  sample_path

# get library size in millions of reads
library_size <- 
  as.integer(sample_table$library_size) %>% 
  set_names(., sample_table$name) %>% 
  divide_by(., 1e6)

# read in the ENSEMBL gene table, take protein coding genes
ensembl_gtf <- 
  read_delim(file = mm10_gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c")) %>% 
  dplyr::filter(str_detect(X9, "protein_coding"))

################################################################################### MAIN CODE
### take transcript with most exons for each gene
# get 5'UTR-s for each transcript
gtf_5utr <- GffToGRanges(ensembl_gtf, "five_prime_utr")

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

# get whole ranges of transcripts for summarasing counts
gtf_range <- unlist(range(gtf_trans))

################################################################################### 
### get count of reads 
# counting overlaps
bamfiles <- Rsamtools::BamFileList(sample_path, yieldSize = 2000000)
se <- GenomicAlignments::summarizeOverlaps(features = gtf_range,
                                           reads = bamfiles,
                                           mode = "Union",
                                           singleEnd = FALSE,
                                           ignore.strand = TRUE)

# get counts
se_counts <- 
  assay(se) %>% 
  as.data.frame() %>% 
  set_colnames(left_join(tibble(ID = str_replace(colnames(.), "_.*", "")), sample_table,  by = "ID") %$% name) %>% 
  tibble::rownames_to_column("feature")

# FPKM normalization
se_fpkm <- 
  cbind(feature = se_counts[, 1], 
        round(data.frame(t(t(se_counts[, -1]) / library_size[colnames(se_counts[, -1])])), 2), 
        feature_width = width(gtf_range) / 1000) %>% 
  mutate_at(.cols = vars(starts_with("s_")), .funs = funs(. / feature_width)) %>% 
  dplyr::select(-feature_width)

# get IDs of top 100 transcripts by expression in 1-cell stage
top_transcriptIDs <- 
  se_fpkm %>% 
  dplyr::select(which(str_detect(colnames(.), "feature|WT"))) %>% 
  mutate(s_1C_WT_mean = rowMeans(dplyr::select(., starts_with("s_")), na.rm = TRUE)) %>% 
  top_n(n = 100, wt = s_1C_WT_mean) %>% 
  dplyr::arrange(desc(s_1C_WT_mean)) %$% 
  feature %>% 
  as.character(.) 

# get ranges of 5'UTRs of top 100 transcripts, write ranges
top_transcript_5utr <-
  gtf_5utr[gtf_5utr$transcript_id %in% top_transcriptIDs] %>%
  as.data.frame() %>% 
  write_delim(path = file.path(outpath, "top100_ensembl_genes_1C_WT_fpkm_expression_5UTR_coordinates.txt"))
  
  
  