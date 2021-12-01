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

library(GenomicAlignments)
library(GenomicRanges)
library(GenomicFeatures)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# path to plasmid coordinates
plasmid_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/Documentation/pCAG-EGFP_MosIR.coordinates.bed"

# mapped path
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Data/Mapped/STAR_mm10.pCAG_EGFP_MosIR.new"

# classified reads path
read_class_path <- list.files(path = file.path(mapped_path, "5_class_reads"), pattern = "read_class.*.txt", full.names = T)

# library size path
library_size_path <- file.path(inpath, "library_size.Eliska_mESC_MosIR.csv")
  
######################################################## READ DATA
# read library size
library_size <- readr::read_csv(library_size_path)

######################################################## MAIN CODE
# set categories order
read_classes <- c("Mos_mRNA", "miRNA.mature.sense", "miRNA.other.sense", "protein_coding.sense", 
                  "rRNA", "SINE", "LINE", "LTR", "other_repeat", "annotated_pseudogene", "other", "not_annotated")

# read class reads table, gets different mismatch counts
read_class_df <-
  purrr::map(read_class_path, readr::read_delim, delim = "\t") %>%
  dplyr::bind_rows(.) %>% 
  dplyr::rename(alignment_count = count, sample_id = sample) %>% 
  dplyr::select(-experiment) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.mis_0\\.18to30nt"))

# get and spread counts 
read_class_counts <- 
  read_class_df %>% 
  tidyr::spread(read_group, alignment_count) %>% 
  dplyr::select(sample_id, read_classes)

# get and spread RPMs
read_class_rpms <- 
  read_class_df %>% 
  dplyr::left_join(., library_size, by = "sample_id") %>% 
  dplyr::mutate(library_size.mis_0.18to30nt = library_size.mis_0.18to30nt / 1e6, 
                rpm = round(alignment_count / library_size.mis_0.18to30nt, 3)) %>% 
  dplyr::select(read_group, rpm, sample_id) %>% 
  tidyr::spread(read_group, rpm) %>% 
  dplyr::select(sample_id, read_classes)

# write to separate sheets in xlsx 
write.xlsx(x = as.data.frame(read_class_counts), 
           file = file.path(outpath, str_c("class_reads.Eliska_mESC_MosIR.hierarchy_counts.xlsx")), 
           sheetName = "allowed_0_mismatches.counts", 
           row.names = FALSE)

# write to separate sheets in xlsx 
write.xlsx(x = as.data.frame(read_class_rpms), 
           file = file.path(outpath, str_c("class_reads.Eliska_mESC_MosIR.hierarchy_counts.xlsx")), 
           sheetName = "allowed_0_mismatches.RPMs", 
           row.names = FALSE, 
           append = TRUE)

