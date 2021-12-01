### INFO: plot coverage over pCMV MosIR EGFP with 21-23 nt reads (normalized)
### DATE: Fri Jul 20 17:34:24 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/siRNA.Valeria/datasets/mouse_spleen.MosIR.small_RNAseq.2021_Jan/Analysis/MosIR_coverage")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(openxlsx)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
### get counts of read alignments with different mismatch cut-off over MosIR in plasmid 
countMosIR <- function(bam_path, which_gr, on_minus_strand = NA){
  
  # get bam name
  bam_name <- basename(bam_path) %>% str_remove(., ".bam")
  
  # read bam
  bam_gr <- 
    GenomicAlignments::readGAlignmentsList(file = bam_path, 
                                           use.names = TRUE, 
                                           param = ScanBamParam(which = which_gr, 
                                                                flag = scanBamFlag(isMinusStrand = on_minus_strand))) %>% 
    unlist(.)
  
  # get counts over MosIR
  mosir_widths <- 
    IRanges::subsetByOverlaps(bam_gr, mosir_gr, ignore.strand = T) %>% 
    width(.)
  
  # return tibble
  count_table <- 
    tibble(alignment_length = mosir_widths) %>% 
    dplyr::group_by(alignment_length) %>% 
    dplyr::summarise(alignment_count = n()) %>% 
    dplyr::ungroup(.) %>% 
    dplyr::right_join(tibble(alignment_length = 12:50), by = "alignment_length") %>% 
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

# MosIR .bed path
bed_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/Documentation/pCAG-EGFP_MosIR.coordinates.bed"

# experiment paths
experiment_path <- "/common/WORK/fhorvat/Projekti/Svoboda/siRNA.Valeria/datasets/mouse_spleen.MosIR.small_RNAseq.2021_Jan"

# mapped path
mapped_path <- file.path(experiment_path, "Data/Mapped/STAR_mm10/5_filtered_reads")

# list bams 
bam_paths <- list.files(path = mapped_path, pattern = ".*\\.bam$", full.names = T)

# library size path
library_size_path <- file.path(mapped_path, "library_sizes.txt")

######################################################## READ DATA
# read bed with coordinates
mosir_gr <- rtracklayer::import.bed(con = bed_path)
mosir_gr <- mosir_gr[mosir_gr$name != "EGFP"]

# read library size df
library_size_df <-  readr::read_delim(file = library_size_path, delim = "\t", col_names = c("sample_id", "library_size"))

######################################################## MAIN CODE
### prepare files
# set whole plasmid coordinates
plasmid_gr <- GenomicRanges::GRanges(seqnames = "pCAG-EGFP_MosIR",
                                     ranges = IRanges(start = 1, end = 6600),
                                     gene_id = "pCAG-EGFP_MosIR")

# clean library size
library_size_tidy <- 
  library_size_df %>% 
  dplyr::filter(str_detect(sample_id, "\\.18to30nt$"))


### get MosIR counts for all samples 
# get counts of MosIR from all bam files in one data.frame
counts_df <- 
  purrr::map(bam_paths, countMosIR, which_gr = plasmid_gr, on_minus_strand = F) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::filter(alignment_length >= 18, alignment_length <= 30) %>% 
  dplyr::mutate(alignment_length = str_c("l_", alignment_length), 
                alignment_length = factor(alignment_length, levels = str_c("l_", 18:30))) %>% 
  dplyr::arrange(alignment_length) %>% 
  dplyr::left_join(., library_size_tidy, by = "sample_id") %>%
  dplyr::mutate(library_size = round((library_size / 1e6), 4),
                rpm = alignment_count / library_size) %>% 
  dplyr::select(sample_id, alignment_length, rpm) %>% 
  tidyr::pivot_wider(., id_cols = sample_id, names_from = alignment_length, values_from = rpm) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.18to30nt$"))
  
# write
wb <- loadWorkbook(file.path(outpath, "MosIR_RPM.18to30nt_reads.summary.xlsx"))
openxlsx::addWorksheet(wb, sheetName = "mis_1")
writeDataTable(wb, sheet = "mis_1", counts_df, colNames = T)
saveWorkbook(wb, file.path(outpath, "MosIR_RPM.18to30nt_reads.summary.xlsx"), overwrite = T)



