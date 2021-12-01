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
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
outpath <- getwd()
library_size_path <- list.files(outpath, "library_sizes.txt", full.names = T)
tracks_path <- list.files(outpath, "log.tracks_URL.txt", full.names = T)
read_class_path <- list.files(outpath, "read_class.*txt", full.names = T)

######################################################## READ DATA

######################################################## MAIN CODE
# get experiment name
experiment_name <- str_remove(outpath, "/Data.*") %>% basename(.)

# library size
library_size <- 
  readr::read_delim(library_size_path, col_names = c("sample_id", "library_size"), delim = "\t") %>% 
  dplyr::mutate(size_select = ifelse(str_detect(sample_id, "21to23nt"), "lib_size.21to23nt", "lib_size.all"), 
                sample_id = str_remove(sample_id, ".21to23nt")) %>% 
  tidyr::spread(size_select, library_size)

# read classes
read_class <- 
  purrr::map(read_class_path, readr::read_delim, delim = "\t") %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::select(-experiment) %>% 
  tidyr::separate()
  dplyr::mutate(size_select = ifelse(str_detect(sample_id, "21to23nt"), "lib_size.21to23nt", "lib_size.all"), 
                sample_id = str_remove_all(sample_id, ".SE|.21to23nt")) %>% 
  tidyr::spread(size_select, library_size) 

# tracks
tracks_tbl <-
  purrr::map(tracks_path, readr::read_delim, col_names = "URL", delim = "\t") %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::filter(str_detect(string = URL, pattern = "http"),
                !str_detect(string = URL, pattern = "bai"), 
                str_detect(string = URL, pattern = str_c(stats_tbl$sample_id, collapse = "|"))) %>%
  dplyr::mutate(sample_id = str_extract(URL, pattern = str_c(stats_tbl$sample_id, collapse = "|")),
                experiment = "ES_DcrTrans_2012",
                scaled = ifelse(str_detect(URL, "scaled"), "RPM_scaled", "raw"),
                size_select = ifelse(str_detect(URL, "21to23nt"), "21to23nt_length", "all_length"), 
                file_type = ifelse(str_detect(basename(URL), "bw"), "coverage", "individual_reads"),
                bw_name = ifelse(scaled == "RPM_scaled", 
                                 yes = str_c(str_remove(sample_id, "s_"), ".bw.scaled"), 
                                 no = str_c( str_remove(sample_id, "s_"), ".bw.raw")), 
                bam_name = str_c(str_remove(sample_id, "s_"), ".bam"), 
                URL = ifelse(file_type == "coverage",
                             yes = str_c("track type=bigWig name=\"", bw_name, "\" bigDataUrl=\"", URL, "\""),
                             no = str_c("track type=bam name=\"", bam_name, "\" bigDataUrl=\"", URL, "\""))) %>%
  dplyr::select(-c(bw_name, bam_name)) %>% 
  tidyr::unite(URL_type, scaled, file_type, size_select, sep = ".") %>%
  tidyr::spread(key = URL_type, value = URL)

# combine and save
stats_and_tracks <-
  right_join(stats_tbl, tracks_tbl, by = "sample_id") %>%
  dplyr::select(experiment, sample_id, lib_size.all, lib_size.21to23nt, contains("all_length"), contains("21to23nt_length")) %T>%
  readr::write_csv(., path = file.path(outpath, str_c("log.", unique(.$experiment), ".stats_and_tracks.csv")))


