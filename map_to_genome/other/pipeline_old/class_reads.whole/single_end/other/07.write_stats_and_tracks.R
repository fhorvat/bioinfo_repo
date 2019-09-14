#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: combine track URLs and read stats table
### DATE: 09. 10. 2017.
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
outpath <- getwd()
stats_path <- file.path(outpath, "log.read_stats.txt")
tracks_path <- file.path(outpath, "log.tracks_URL.txt")

######################################################## READ DATA

######################################################## MAIN CODE
# stats
stats_tbl <- readr::read_delim(stats_path, delim = "\t")

# tracks
tracks_tbl <-
  readr::read_delim(file = tracks_path, col_names = F, delim = "\t") %>%
  dplyr::select(URL = X1) %>%
  dplyr::filter(str_detect(string = URL, pattern = "http"),
                !str_detect(string = URL, pattern = "bai"), 
                str_detect(string = URL, pattern = str_c(stats_tbl$sample_id, collapse = "|"))) %>%
  dplyr::mutate(sample_id = str_extract(URL, pattern = str_c(stats_tbl$sample_id, collapse = "|")),
                experiment = str_replace(URL, pattern = str_c(str_c("/", unique(sample_id), ".*"), collapse = "|"), replacement = "") %>% basename(.),
                scaled = ifelse(str_detect(URL, "scaled"), "RPM_scaled", "raw"),
                file_type = ifelse(str_detect(basename(URL), "bw"), "coverage", "individual_reads"),
                bw_name = ifelse(test = (scaled == "RPM_scaled"), 
                                 yes = str_c(experiment, ".", str_replace(sample_id, "s_", ""), ".scaled"), 
                                 no = str_c(experiment, ".", str_replace(sample_id, "s_", ""), ".raw")), 
                bam_name = str_c(experiment, ".", str_replace(sample_id, "s_", "")), 
                URL = ifelse(file_type == "coverage",
                             yes = str_c("track type=bigWig name=\"", bw_name, "\" bigDataUrl=\"", URL, "\""),
                             no = str_c("track type=bam name=\"", bam_name, "\" bigDataUrl=\"", URL, "\""))) %>%
  dplyr::select(-c(bw_name, bam_name)) %>% 
  tidyr::unite(scale_type, scaled, file_type, sep = ".") %>%
  tidyr::spread(key = scale_type, value = URL) %>%
  dplyr::select(sample_id, experiment, raw.coverage, RPM_scaled.coverage, raw.individual_reads)

# combine and save
stats_and_tracks <-
  right_join(stats_tbl, tracks_tbl, by = "sample_id") %>%
  dplyr::select(experiment, everything()) %T>%
  readr::write_csv(., path = file.path(outpath, str_c("log.", unique(.$experiment), ".stats_and_tracks.csv")))


