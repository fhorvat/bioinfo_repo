#!/home/students/fhorvat/R/bin/Rscript
### INFO: checks genomic position for overlap with annotated and non-canonical splice junctions from STAR output
### DATE: 14. 09. 2017.
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Split_RNAseq_2017/mouse/Data/Mapped/STAR_mm10/SINE_region")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(tibble)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)
library(BiocParallel)

######################################################## PATH VARIABLES
repeatmasker_path <- "/common/WORK/fhorvat/reference/mouse/mm10/UCSC/UCSC_repeatMasker_withLTRsubclass.txt.gz"
outpath <- getwd()
blast_list <- list.files(path = outpath, pattern = "*blastn.tbl$", full.names = T)
file <- blast_list[1]
file_name <- 
  basename(file) %>% 
  str_replace(., pattern = ".blastn.tbl", "")

# repeatMasker
rmsk <-
  readr::read_delim(file = repeatmasker_path, delim = "\t", col_names = T) 
  # dplyr::mutate(full_pos = str_c(seqnames, ":", start, "-", end, "|", strand, "|", repClass, "|", repName)) %>%
  # GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "GffToGRanges.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))

######################################################## FUNCTIONS


######################################################## READ DATA
# read BLAST format 6 table
blast6 <- 
  readr::read_delim(file = file, delim = "\t",
                    col_names = c("query_id", "subject_id", "identity_perc", "alignment_length", 
                                  "mismatches", "gap_open", 
                                  "query_start", "query_end", "subject_start", "subject_end", 
                                  "e_value", "bit_score"))

######################################################## MAIN CODE
# filter gap in reads
blast_gap_reads <- 
  blast6 %>% 
  dplyr::filter(subject_id == "mm10_dna") %$%
  query_id

# get whole alignments
blast_whole <- 
  blast6 %>% 
  dplyr::filter(subject_id != "mm10_dna", 
                query_id %in% blast_gap_reads, 
                alignment_length > 44) %$%
  query_id %>%
  readr::write_lines(x = ., path = file.path(outpath, stringr::str_c(file_name, ".readsIDs.txt")))

# check for repeat chr19 5801969 5802116	-	B1_Mus1	SINE Alu
rmsk_filt <-
  rmsk %>% 
  dplyr::filter(repClass == "B1_Mus1") %>% 
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)
  