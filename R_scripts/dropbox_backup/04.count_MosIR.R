### INFO: plot coverage over pCMV MosIR EGFP with 21-23 nt reads (normalized)
### DATE: Fri Jul 20 17:34:24 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/MosIR_expression")

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
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

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
mosir_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/Documentation/pCAG-EGFP_MosIR.coordinates.bed"

######################################################## READ DATA
# read bed with coordinates
mosir_gr <- rtracklayer::import.bed(con = mosir_path)
mosir_gr <- mosir_gr[mosir_gr$name != "EGFP"]

######################################################## MAIN CODE
# set whole plasmid coordinates
plasmid_gr <- GenomicRanges::GRanges(seqnames = "pCAG-EGFP_MosIR", 
                                     ranges = IRanges(start = 1, end = 10000), 
                                     gene_id = "pCAG-EGFP_MosIR")

### get MosIR counts for all samples in each experiment 
# experiment paths
experiment_paths <- 
  c("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/SOLiD/T3T_DcrTrans_2011/Data/Mapped/Shrimp_mm10", 
    "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/SOLiD/ES_DcrTrans_2012/Data/Mapped/Shrimp_mm10", 
    "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/T3T_mESC_MosIR.2016/Data/Mapped/STAR_mm10",
    "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Data/Mapped/STAR_mm10.pCAG_EGFP_MosIR") %>% 
  magrittr::set_names(., str_remove(., "Data.*$") %>% basename(.))

# list of experiments
experiment_list <- names(experiment_paths)

# loop through experiments
for(experiment in experiment_list){
  
  # get mapped path 
  mapped_path <- experiment_paths[experiment]
  
  # list bams 
  bam_paths <- list.files(path = mapped_path, pattern = ".*mis_.*bam$", full.names = T)
  
  # get counts of MosIR from all bam files in one data.frame, sum 0 mismatches with 1 and 2 mismatches counts
  counts_df <- 
    purrr::map(bam_paths, countMosIR, which_gr = plasmid_gr, on_minus_strand = F) %>% 
    dplyr::bind_rows(.) %>% 
    tidyr::separate(sample_id, c("sample_id", "mismatch"), sep = ".mis_") %>% 
    dplyr::mutate(mismatch = str_c("mis_", mismatch)) %>% 
    tidyr::spread(mismatch, alignment_count) %>% 
    dplyr::mutate_at(vars(starts_with("mis")), funs(replace(., is.na(.), 0))) %>% 
    dplyr::mutate(mis_1 = mis_1 + mis_0, 
                  mis_2 = mis_2 + mis_1) %>% 
    dplyr::arrange(sample_id, alignment_length) %>% 
    dplyr::filter(alignment_length >= 18, alignment_length <= 30) %>% 
    tidyr::gather(mismatch_allowed, alignment_count, mis_0:mis_2) %>% 
    tidyr::unite(temp, mismatch_allowed, alignment_length) %>%
    tidyr::spread(temp, alignment_count) %>%
    dplyr::mutate_at(vars(starts_with("mis_")), funs(replace(., is.na(.), 0)))
  
  # write to separate sheets in xlsx 
  write.xlsx(x = counts_df %>% 
               dplyr::select_at(vars("sample_id", starts_with("mis_0"))) %>% 
               magrittr::set_colnames(., str_replace(colnames(.), "mis_.*_", "l_")) %>% 
               as.data.frame(.), 
             file = file.path(outpath, str_c("mosIR.", experiment, ".counts_summary.xlsx")), 
             sheetName = "allowed_0_mismatches", 
             row.names = FALSE)
  
  write.xlsx(x = counts_df %>% 
               dplyr::select_at(vars("sample_id", starts_with("mis_1"))) %>% 
               magrittr::set_colnames(., str_replace(colnames(.), "mis_.*_", "l_")) %>% 
               as.data.frame(.), 
             file = file.path(outpath, str_c("mosIR.", experiment, ".counts_summary.xlsx")), 
             sheetName = "allowed_1_mismatches", 
             append = TRUE, 
             row.names = FALSE)
  
  write.xlsx(x = counts_df %>% 
               dplyr::select_at(vars("sample_id", starts_with("mis_2"))) %>% 
               magrittr::set_colnames(., str_replace(colnames(.), "mis_.*_", "l_")) %>% 
               as.data.frame(.), 
             file = file.path(outpath, str_c("mosIR.", experiment, ".counts_summary.xlsx")), 
             sheetName = "allowed_2_mismatches", 
             append = TRUE, 
             row.names = FALSE)
  
}
  
