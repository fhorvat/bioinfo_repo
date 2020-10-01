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
# set outpath, get paths for gtf and bam files
outpath <- getwd()

# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# arguments from command line
experiment <- args$experiment
experiment_name <- args$experiment_name

# stats and tracks paths
stats_path <- file.path(outpath, "bismark_summary_report.txt")
tracks_path <- file.path(outpath, "log.tracks_URL.txt")

######################################################## READ DATA

######################################################## MAIN CODE
# stats
stats_tbl <- 
  readr::read_delim(stats_path, delim = "\t") %>% 
  dplyr::select(-c(`No Genomic Sequence`, `Duplicate Reads (removed)`, `Unique Reads (remaining)`)) %>% 
  dplyr::rename(sample_id = File) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "_bismark_bt2\\.bam|_bismark_bt2_pe\\.bam"))

# tracks
tracks_tbl <-
  readr::read_delim(file = tracks_path, col_names = F, delim = "\t") %>%
  dplyr::select(URL = X1) %>%
  dplyr::filter(str_detect(string = URL, pattern = "\\.bw$|\\.bam$"),
                str_detect(string = URL, pattern = "http"),
                !str_detect(string = URL, pattern = "bai"),
                str_detect(string = URL, pattern = str_c(stats_tbl$sample_id, collapse = "|"))) %>%
  dplyr::mutate(sample_id = basename(URL) %>% str_remove_all(., "\\_bismark_bt2\\.bam$|\\.bw$|\\.scaled|_bismark_bt2_pe\\.bam|_bismark_bt2_pe"),
                experiment = experiment,
                file_type = ifelse(test = str_detect(basename(URL), "bw"),
                                   yes = "methylation_coverage",
                                   no = "individual_reads"),
                bw_name = str_remove(sample_id, "^s_"),
                bam_name = str_c(str_remove(sample_id, "^s_"), "bam", sep = "."),
                URL = ifelse(test = (file_type == "methylation_coverage"),
                             yes = str_c("track type=bigWig name=\"", bw_name, "\" bigDataUrl=\"", URL, "\""),
                             no = str_c("track type=bam name=\"", bam_name, "\" bigDataUrl=\"", URL, "\""))) %>%
  dplyr::select(-c(bw_name, bam_name))

# reports
reports_tbl <- 
  readr::read_delim(file = tracks_path, col_names = F, delim = "\t") %>%
  dplyr::select(report_URL = X1) %>%
  dplyr::filter(str_detect(string = report_URL, pattern = "\\_report.html$"),
                str_detect(string = report_URL, pattern = str_c(stats_tbl$sample_id, collapse = "|"))) %>% 
  dplyr::mutate(sample_id = basename(report_URL) %>% str_remove_all(., "_bismark_bt2_.*_report\\.html"))

# full report
report_all_tbl <- 
  readr::read_delim(file = tracks_path, col_names = F, delim = "\t") %>% 
  dplyr::select(sample_id = X1) %>% 
  dplyr::filter(str_detect(sample_id, "bismark_summary_report.html"))

# create table with all possible columns
tracks_placeholder <-
  tibble(sample_id = rep(unique(tracks_tbl$sample_id), each = 2),
         file_type = rep(c("methylation_coverage", "individual_reads"), length(unique(tracks_tbl$sample_id))))

# join placeholder with table
tracks_tbl_tidy <-
  tracks_tbl %>%
  dplyr::select(-experiment) %>%
  dplyr::right_join(., tracks_placeholder, by = c("sample_id", "file_type")) %>%
  tidyr::spread(key = file_type, value = URL) %>%
  dplyr::mutate(experiment = experiment)

# combine and save
stats_and_tracks <-
  right_join(stats_tbl, tracks_tbl_tidy, by = "sample_id") %>%
  right_join(., reports_tbl, by = "sample_id") %>% 
  dplyr::select(experiment, sample_id, methylation_coverage, individual_reads, report_URL, everything()) %>% 
  bind_rows(., report_all_tbl) %T>%
  readr::write_csv(., path = file.path(outpath, str_c("log.", experiment, ".stats_and_tracks.csv")))
