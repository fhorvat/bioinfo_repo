### INFO: R Script
### DATE: 22. 06. 2017. 
### AUTHOR: Filip Horvat
### PATH: 
rm(list = ls()); gc()
# options(bitmapType = 'cairo')

################################################################################### WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/SMRT_oocytes_2017/Analysis/polyA_length")

################################################################################### LIBRARIES
# data shaping
library(dplyr)
library(data.table)
library(stringr)
library(readr)
library(magrittr)
library(ggplot2)
library(reshape2)
library(tibble)

# genomics
library(ShortRead)
library(org.Mm.eg.db)
library(EnsDb.Mmusculus.v79)

################################################################################### PATH VARIABLES
fastq_path <- "/common/WORK/fhorvat/Projekti/Svoboda/SMRT_oocytes_2017/Data/Raw/Cleaned/02_demultiplexed/bbtools_seal"
fastq_files <- list.files(fastq_path, pattern = "[B6|DBA]_[GV|MII].*fastq$", full.names = T)

bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/SMRT_oocytes_2017/Data/Mapped/STARlong_mm10"
bam_files <- list.files(bam_path, pattern = "[B6|DBA]_[GV|MII].*bam$", full.names = T)
  
outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/SMRT_oocytes_2017/Analysis/polyA_length"
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
mm10_gtf_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Ensembl_GRCm38.86.20161128.gtf.gz"

################################################################################### SOURCE FILES
source(file.path(lib_path, "GffToGRanges.R"))

################################################################################### FUNCTIONS

################################################################################### SCRIPT PARAMS

################################################################################### TABLES 
# read in the ENSEMBL gene table, take protein coding genes
ensembl_gtf <- read_delim(file = mm10_gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c")) 

# genes of interest - Mos, cyclin B1, dcp1a
genes_list <- 
  tibble::tibble(gene_symbol = c("Mos", "Ccnb1", "Dcp1a")) %>% 
  dplyr::mutate(ensembl_txid = mapIds(EnsDb.Mmusculus.v79, keys = gene_symbol, column = "TXID", keytype = "SYMBOL", multiVals = "first"), 
                ensembl_geneid = mapIds(EnsDb.Mmusculus.v79, keys = gene_symbol, column = "GENEID", keytype = "SYMBOL", multiVals = "first"))

# read mapped .bam files (GV and MII)
bam_list <- 
  lapply(bam_files, GenomicAlignments::readGAlignments, use.names = T) %>% 
  set_names(str_replace_all(string = bam_files, pattern = "\\/.*\\/|_Aligned.sortedByCoord.out.bam", replacement = ""))

################################################################################### MAIN CODE
# get exons for each transcript
gtf_trans <- GffToGRanges(ensembl_gtf, "exon")
gtf_trans <- split(gtf_trans, gtf_trans$transcript_id)
gtf_trans_filter <- gtf_trans[names(gtf_trans) %in% genes_list$ensembl_txid]

### loop through fastq files, get length of polyA tails in reads which are mapped to reads of interest
all_mapped_reads <- lapply(X = 1:length(bam_list), FUN = function(X){
  
  # find overlaps
  overlaps_exons <- GenomicRanges::findOverlaps(query = gtf_trans_filter, subject = bam_list[[X]], ignore.strand = T) 
  
  # get read ID's
  mapped_reads <- 
    tibble(ensembl_txid = names(gtf_trans_filter[queryHits(overlaps_exons)]), 
           read_id = names(bam_list[[X]][subjectHits(overlaps_exons)]), 
           bam_file = names(bam_list)[X]) %>% 
    dplyr::left_join(genes_list %>% dplyr::select(gene_symbol, ensembl_txid), .) %>% 
    dplyr::select(read_id, gene_symbol, bam_file) 
  
  # read fastq file, get sequences and read ID's
  fastq_all <- ShortRead::readFastq(dirPath = fastq_files[X])
  fastq_seq <- ShortRead::sread(fastq_all)
  fastq_id <- ShortRead::id(fastq_all)
  
  # filter reads mapped to genes of interests from fastq, convert to character vectors
  fastq_seq_filter <- sapply(fastq_seq[which(fastq_id %in% mapped_reads$read_id)], toString)
  
  # get length of shortest polyA tail in sequences
  polyA_min_length <-
    stringr::str_extract_all(string = fastq_seq_filter, pattern = "A{5,}|T{5,}") %>%
    sapply(., function(x) max(nchar(x))) %>%
    min
  
  # get polyA from each sequence (with 10% edit distance), get length
  polyA_length <-
    regmatches(fastq_seq_filter, aregexec(pattern = str_c("A{", polyA_min_length, ",}|T{", polyA_min_length, ",}"),
                                          text = fastq_seq_filter,
                                          max.distance = 0.1)) %>%
    lapply(., nchar) %>%
    lapply(., function(x) replace(x, length(x) == 0, 0)) %>%
    unlist(.)
  
  # add length to genes list
  mapped_reads %<>%
    mutate(polyA_length = polyA_length) %>% 
    dplyr::select(-read_id)
  
  return(mapped_reads)
  
}) %>% 
  dplyr::bind_rows(.) 

### plot as jitterplot + boxplot
ggplot(data = all_mapped_reads, aes(bam_file, polyA_length)) +
  geom_boxplot(show.legend = F, outlier.shape = NA) +
  geom_jitter(aes(colour = bam_file), size = 1.5, width = 0.25, show.legend = F) +
  facet_grid(gene_symbol ~ .) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("sample") +
  ylab("polyA length") +
  ggsave(filename = file.path(outpath, "polyA_length_3genes_MII_GV_revComp_reads.pdf"), width = 25, height = 15)
  