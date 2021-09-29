#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: combine track URLs and read stats table
### DATE: 09. 10. 2017.
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd(".")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set outpath
outpath <- getwd()

# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# arguments from command line
experiment <- args$experiment
experiment_name <- args$experiment_name
log_path <- args$log_path

# tracks paths
tracks_path <- file.path(outpath, "log.tracks_URL.txt")

######################################################## READ DATA

######################################################## MAIN CODE
# stats
stats_tbl <- 
  readr::read_delim(log_path, delim = "\t") %>% 
  dplyr::mutate(subset = str_extract(sample_id, "(?<=\\.[S,P]E\\.).*$"),
                subset = ifelse(is.na(subset), "whole_set", subset),
                sample_id = str_remove(sample_id, str_c("\\.", subset)) %>% str_remove(., "_r[0-9]+(?=\\.[S,P]E$)")) %>% 
  dplyr::group_by(sample_id, subset) %>% 
  dplyr::summarise_all(sum) %>% 
  dplyr::ungroup(.) %>% 
  tidyr::unite(col = "sample_id", sample_id, subset, sep = ".") %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.whole_set$"))

# tracks
tracks_tbl <-
  readr::read_delim(file = tracks_path, col_names = F, delim = "\t") %>%
  dplyr::select(URL = X1) %>%
  dplyr::filter(str_detect(string = URL, pattern = "\\.bw$|\\.bam$"),
                str_detect(string = URL, pattern = "http"),
                !str_detect(string = URL, pattern = "bai")) %>% 
  dplyr::filter(str_detect(string = URL, pattern = str_c(stats_tbl$sample_id, collapse = "|"))) %>%
  dplyr::mutate(sample_id = basename(URL) %>% str_remove_all(., "\\.bam$|\\.bw$|\\.scaled|\\_19to32"),
                experiment = experiment,
                experiment_short = experiment_name,
                scaled = ifelse(str_detect(URL, "scaled"), "RPM_scaled", "raw"),
                scaled = ifelse(str_detect(URL, "_19to32"), str_c(scaled, "_19to32"), scaled),
                file_type = ifelse(test = str_detect(basename(URL), "bw"),
                                   yes = "coverage",
                                   no = "individual_reads"),
                bw_name = str_c(str_remove_all(sample_id, "^s_|\\.SE|''.PE"), scaled, sep = "."),
                bam_name = str_c(str_remove_all(sample_id, "^s_|\\.SE|\\.PE"), sep = "."),
                URL = ifelse(test = (file_type == "coverage"),
                             yes = str_c("track type=bigWig name=\"", bw_name, "\" bigDataUrl=\"", URL, "\""),
                             no = str_c("track type=bam name=\"", bam_name, "\" bigDataUrl=\"", URL, "\""))) %>%
  dplyr::select(-c(bw_name, bam_name)) %>%
  tidyr::unite(scale_type, scaled, file_type, sep = ".") %>%
  tidyr::spread(key = scale_type, value = URL) %>%
  dplyr::select(sample_id, experiment, contains("coverage"), raw.individual_reads)

# combine and save
stats_and_tracks <-
  right_join(stats_tbl, tracks_tbl, by = "sample_id") %>%
  dplyr::mutate(subset = str_extract(sample_id, "(?<=\\.SE\\.).*$"),
                subset = ifelse(is.na(subset), "whole_set", subset),
                sample_id = str_remove(sample_id, str_c("\\.", subset))) %>%
  dplyr::select(experiment, sample_id, subset, contains("coverage"), raw.individual_reads, everything()) %T>%
  readr::write_csv(., file = file.path(outpath, str_c("log", experiment, "stats_and_tracks.csv", sep = ".")))
