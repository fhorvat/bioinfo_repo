### INFO: plot coverage over pCMV MosIR EGFP with 21-23 nt reads (normalized)
### DATE: Wed Nov 28 15:03:32 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Analysis")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(xlsx)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
### get counts of read alignments with different mismatch cut-off over MosIR in plasmid 
countIR <- function(bam_path){
  
  # get bam name
  bam_name <- basename(bam_path) %>% str_remove(., "\\.SE\\.mis_0\\.18to30nt\\.bam")

  # read bam
  bam_gr <- 
    GenomicAlignments::readGAlignmentsList(file = bam_path, 
                                           use.names = TRUE,
                                           param = ScanBamParam(which = ir_gr, 
                                                                flag = scanBamFlag(isMinusStrand = NA))) %>% 
    unlist(.)
  
  # return tibble
  count_table <- 
    tibble(read_id = names(bam_gr), 
           alignment_length = width(bam_gr)) %>% 
    dplyr::group_by(read_id) %>% 
    dplyr::summarise(alignment_length = max(alignment_length)) %>% 
    dplyr::group_by(alignment_length) %>% 
    dplyr::summarise(alignment_count = n()) %>% 
    dplyr::ungroup(.) %>% 
    dplyr::right_join(tibble(alignment_length = rep(18:30)), by = "alignment_length") %>% 
    dplyr::mutate(alignment_count = replace(alignment_count, is.na(alignment_count), 0)) %>% 
    dplyr::mutate(sample_id = bam_name)
  
  # return
  return(count_table)
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# get path of inverted repeat coordinates from mRNA
ir_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Analysis/MosIR.coordinates.csv"

# get mapped path 
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Data/Mapped/STAR_mm10.pCAG_EGFP_MosIR.new"

# list bams 
bam_paths <- list.files(path = mapped_path, pattern = ".*bam$", full.names = T)

# library size path
library_size_path <- file.path(inpath, "library_size.Eliska_mESC_MosIR.csv")

######################################################## READ DATA
# read coordinates of inverted repeats from mRNA
ir_gr <- readr::read_csv(file = ir_path)

# read library size
library_size <- readr::read_csv(library_size_path)

######################################################## MAIN CODE
# prepare mRNA IR coordinates on plasmids
ir_gr %<>% GenomicRanges::GRanges(.)

# get counts of MosIR from all bam files in one data.frame
counts_df <- 
  purrr::map(bam_paths, countIR) %>% 
  dplyr::bind_rows(.) 

# spread counts
mosir_counts <- 
  counts_df %>% 
  tidyr::spread(alignment_length, alignment_count) %>% 
  magrittr::set_colnames(., c("sample_id", str_c("l_", colnames(.)[2:ncol(.)])))

# get and spread RPMs
mosir_rpms <- 
  counts_df %>% 
  dplyr::left_join(., library_size, by = "sample_id") %>% 
  dplyr::mutate(library_size.mis_0.18to30nt = library_size.mis_0.18to30nt / 1e6, 
                rpm = round(alignment_count / library_size.mis_0.18to30nt, 3)) %>% 
  dplyr::select(alignment_length, rpm, sample_id) %>% 
  tidyr::spread(alignment_length, rpm) %>% 
  magrittr::set_colnames(., c("sample_id", str_c("l_", colnames(.)[2:ncol(.)])))

### write to separate sheets in xlsx 
# counts
write.xlsx(x = as.data.frame(mosir_counts), 
           file = file.path(outpath, "mosIR.Eliska_mESC_MosIR.read_length.xlsx"), 
           sheetName = "allowed_0_mismatches.counts", 
           row.names = FALSE)

# RPM
write.xlsx(x = as.data.frame(mosir_rpms), 
           file = file.path(outpath, "mosIR.Eliska_mESC_MosIR.read_length.xlsx"), 
           sheetName = "allowed_0_mismatches.RPMs", 
           row.names = FALSE, 
           append = TRUE)

