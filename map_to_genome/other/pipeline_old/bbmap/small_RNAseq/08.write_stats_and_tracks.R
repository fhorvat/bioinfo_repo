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
library(ggplot2)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
outpath <- getwd()
stats_path <- list.files(file.path(outpath, "5_class_reads"), pattern = "read_class\\..*\\.txt", full.names = T)
tracks_path <- file.path(outpath, "log.tracks_URL.txt")

######################################################## READ DATA

######################################################## MAIN CODE
# set experiment and experiment name
experiment <- "DicerX.spleen"
experiment_name <- "2019_Nov"

# set class hierarchy
class_hier <- c("miRNA.mature.sense", "miRNA.other.sense", "protein_coding.sense", "rRNA", "SINE", "LINE", "LTR", "other_repeat", "annotated_pseudogene", "other", "not_annotated")

# stats
stats_tbl <- purrr::map(stats_path, function(path){
  
  # read data
  suppressMessages(readr::read_delim(file = path, delim = "\t"))
  
}) %>%
  bind_rows(.) %>%
  dplyr::select(-experiment) %>%
  dplyr::mutate(read_group = factor(read_group, levels = class_hier)) %>%
  tidyr::spread(read_group, count) %>%
  dplyr::rename(sample_id = sample)

# tracks
tracks_tbl <-
  readr::read_delim(file = tracks_path, col_names = F, delim = "\t") %>%
  dplyr::select(URL = X1) %>%
  dplyr::filter(str_detect(string = URL, pattern = "\\.bw$|\\.bam$"),
                str_detect(string = URL, pattern = "http"),
                !str_detect(string = URL, pattern = "bai")) %>%
  dplyr::mutate(sample_id = basename(URL) %>% str_remove_all(., "\\.bam$|\\.bw$|\\.scaled"),
                experiment = experiment,
                experiment_short = experiment_name,
                scaled = ifelse(str_detect(URL, "scaled"), "RPM_scaled", "raw"),
                file_type = ifelse(test = str_detect(basename(URL), "bw"),
                                   yes = "coverage",
                                   no = "individual_reads"),
                bw_name = ifelse(test = (scaled == "RPM_scaled"),
                                 yes = str_c(str_remove_all(sample_id, "^s_|\\.SE|\\.PE"), experiment_short, "scaled", sep = "."),
                                 no = str_c(str_remove_all(sample_id, "^s_|\\.SE|''.PE"), experiment_short, "raw", sep = ".")),
                bam_name = str_c(str_remove_all(sample_id, "^s_|\\.SE|\\.PE"), experiment_short, sep = "."),
                URL = ifelse(test = (file_type == "coverage"),
                             yes = str_c("track type=bigWig name=\"", bw_name, "\" bigDataUrl=\"", URL, "\""),
                             no = str_c("track type=bam name=\"", bam_name, "\" bigDataUrl=\"", URL, "\""))) %>%
  dplyr::select(-c(bw_name, bam_name)) %>%
  tidyr::unite(scale_type, scaled, file_type, sep = ".") %>%
  tidyr::spread(key = scale_type, value = URL) %>%
  dplyr::select(sample_id, experiment, contains("coverage"), raw.individual_reads)

# combine and save
stats_and_tracks <-
  inner_join(stats_tbl, tracks_tbl, by = "sample_id") %>%
  dplyr::mutate(subset = str_extract(sample_id, "(?<=\\.SE\\.).*$"),
                subset = ifelse(is.na(subset), "whole_set", subset),
                sample_id = str_remove(sample_id, str_c("\\.", subset))) %>%
  dplyr::select(experiment, sample_id, subset, contains("coverage"), raw.individual_reads, everything()) %>% 
  dplyr::select(-contains("raw.coverage")) %T>%
  readr::write_csv(., path = file.path(outpath, str_c("log.", experiment, ".stats_and_tracks.csv")))

