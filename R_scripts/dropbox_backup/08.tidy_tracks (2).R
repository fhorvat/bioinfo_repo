#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: tidy bw and bam tracks URLs and save
### DATE: Tue Feb 19 16:46:24 2019
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("%OUT_PATH")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# output path
outpath <- getwd()

# tracks path
tracks_path <- file.path(outpath, "log.tracks_URL.txt")

######################################################## READ DATA

######################################################## MAIN CODE
# tidy and save tracks
tracks_tbl <-
  readr::read_delim(file = tracks_path, col_names = F, delim = "\t") %>%
  dplyr::select(URL = X1) %>%
  dplyr::filter(str_detect(string = URL, pattern = "\\.bw$|\\.bam$"), 
                str_detect(string = URL, pattern = "http"),
                !str_detect(string = URL, pattern = "bai")) %>%
  dplyr::mutate(sample_id = basename(URL) %>% str_remove_all(., "\\.bam$|\\.bw$"),
                experiment = dirname(URL) %>% basename(.),
                mismatch_type = ifelse(str_detect(sample_id, "perfect"), "perfect", "mismatches"),
                file_type = ifelse(str_detect(basename(URL), "bw"), "coverage", "individual_reads"),
                scale_type = ifelse(str_detect(sample_id, "scaled"), "scaled", "raw"),
                bw_name = str_c(experiment, str_remove(sample_id, "^s_"), "bw", sep = "."),
                bam_name = str_c(experiment, str_remove(sample_id, "^s_"), "bam", sep = "."),
                URL = ifelse(file_type == "coverage",
                             yes = str_c("track type=bigWig name=\"", bw_name, "\" bigDataUrl=\"", URL, "\""),
                             no = str_c("track type=bam name=\"", bam_name, "\" bigDataUrl=\"", URL, "\"")),
                sample_id = str_remove_all(sample_id, ".perfect|.scaled")) %>%
  dplyr::select(-c(bw_name, bam_name)) %>%
  tidyr::unite(mismatch_scale_file_type, mismatch_type, scale_type, file_type, sep = ".") %>%
  tidyr::spread(key = mismatch_scale_file_type, value = URL) %>%
  dplyr::select(sample_id, experiment, mismatches.raw.coverage, mismatches.scaled.coverage, mismatches.raw.individual_reads, 
                perfect.raw.coverage, perfect.scaled.coverage, perfect.raw.individual_reads) %T>%
  readr::write_csv(., path = file.path(outpath, str_c("log.", unique(.$experiment), ".tracks.csv")))
