#!/home/students/fhorvat/R/bin/Rscript
### INFO: R Script
### DATE: 20. 4. 2017.  
### AUTHOR: Filip Horvat
### PATH: 
### DECRIPTION: reads .fastq files and splits them by list of barcodes
rm(list = ls()); gc()

################################################################################### WORKING DIRECTORY
################################################################################### 
setwd("/common/WORK/fhorvat/Projekti/Svoboda/SMRT_oocytes_2017/Data/Raw/Cleaned/R_demultiplexed")

################################################################################### LIBRARIES
################################################################################### 
library(dplyr)
library(stringr)
library(readr)
library(magrittr)
library(tibble)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(seqinr)
library(Biostrings)
library(ShortRead)

################################################################################### PATH VARIABLES
################################################################################### 
inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/SMRT_oocytes_2017/Data/Raw/Links"

outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/SMRT_oocytes_2017/Data/Raw/Cleaned/R_demultiplexed/output"

barcode_path <- "/common/WORK/fhorvat/Projekti/Svoboda/SMRT_oocytes_2017/Data/Raw/documentation/used_barcodes.fasta"
  
fastq_path <- list.files(path = inpath, pattern = ".fastq", full.names = T)

################################################################################### SOURCE FILES
################################################################################### 

################################################################################### FUNCTIONS
################################################################################### 

################################################################################### SCRIPT PARAMS
################################################################################### 

################################################################################### TABLES
################################################################################### 

################################################################################### MAIN CODE
################################################################################### 
# get barcode sequences
barcode_seq <- readDNAStringSet(filepath = barcode_path)
barcode_seq_df <- tibble(barcode = as.character(barcode_seq), 
                         sample = names(barcode_seq))

# read .fastq, split by barcodes, write to splitted .fastq files
invisible(lapply(fastq_path, function(x){
  
  # get sequences 
  reads_fastq <- ShortRead::readFastq(dirPath = x)
  reads_seq <- ShortRead::sread(reads_fastq)
  reads_id <- as.character(ShortRead::id(reads_fastq))
  
  # get data frame with read id and matching barcode
  barcoded_reads <- lapply(barcode_seq, function(x){
    
    # get info about barcode pattern in begining/end
    barcode_forward <- Biostrings::vmatchPattern(pattern = x, subject = reads_seq, max.mismatch = 2)
    barcode_reverse <- Biostrings::vmatchPattern(pattern = reverseComplement(x), subject = reads_seq, max.mismatch = 2)
    
    # get start/end indices (replace NULL with NA so unlist keeps emtpy indices, if there is more than one match keep only first one)
    start_forward <- startIndex(barcode_forward)
    start_forward[sapply(start_forward, is.null)] <- NA
    start_forward <- lapply(start_forward, function(x) x[1])
    start_forward <- unlist(start_forward)
    
    start_reverse <- endIndex(barcode_reverse) 
    start_reverse[sapply(start_reverse, is.null)] <- NA
    start_reverse <- lapply(start_reverse, function(x) x[length(x)])
    start_reverse <- unlist(start_reverse)
    
    # get table with all info, filter
    barcode_distribution <- 
      tibble(match_forward = elementNROWS(barcode_forward), 
             match_reverse = elementNROWS(barcode_reverse), 
             start_forward = start_forward, 
             start_reverse = start_reverse, 
             read_length = width(reads_seq), 
             read_id = reads_id, 
             barcode = as.character(x)) %>%
      dplyr::filter(!(match_forward == 1 & match_reverse == 1)) %>% 
      dplyr::filter(match_forward == 1 | match_reverse == 1) %>% 
      dplyr::filter((match_forward == 1 & start_forward < 100) | (match_reverse == 1 & (read_length - start_reverse) < 100)) %>% 
      dplyr::select(read_id, barcode)
    
    return(barcode_distribution)
    
  }) %>% 
    dplyr::bind_rows(.) %>% 
    dplyr::group_by(read_id) %>% 
    dplyr::mutate(count = n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(count == 1) %>% 
    dplyr::select(-count) %>% 
    dplyr::left_join(barcode_seq_df, by = "barcode")
  
  # write statistics
  readr::write_lines(x = str_c(basename(x), sum(reads_id %in% barcoded_reads$read_id), sum(!(reads_id %in% barcoded_reads$read_id)), sep = "\t"),
                     path = file.path(outpath, stringr::str_c(basename(x), "_statistics.txt", sep = "")))
  
  # write unmatched reads to unmatched.fastq
  ShortRead::writeFastq(object = reads_fastq[!(reads_id %in% barcoded_reads$read_id)], 
                        file = file.path(outpath, "unmatched.fastq"), 
                        compress = F, 
                        mode = "a")
  
  # split by barcode
  barcoded_reads_list <- plyr::dlply(barcoded_reads, "barcode", identity)
  
  # write as .fastq with append mode
  invisible(lapply(barcoded_reads_list, function(x){
    
    # get sample name
    sample_name <- x$sample[1]
    
    # write matched reads as .fastq
    ShortRead::writeFastq(object = reads_fastq[reads_id %in% x$read_id], 
                          file = file.path(outpath, stringr::str_c(sample_name, ".fastq", sep = "")), 
                          compress = F, 
                          mode = "a")
    
  }))
  
  # print message
  cat("done writing", basename(x), "\n")
  
}))




